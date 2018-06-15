from __future__ import print_function
from argparse import ArgumentParser
import os
import sys
import sqlite3
from contextlib import closing

# keys for db interface
# RUN_NAME = 'run_id'
# FAST5_LOCATION = 'fast5_location'
# READ_ID = 'read_id'
# FAST5_ROOT = 'fast5_root'

RUN_NAME = 'run_id'
FAST5_LOCATION = 'file_location'
READ_ID = 'read_id'
FAST5_ROOT = 'fast5_root'

def parse_args():
    parser = ArgumentParser(description="Build indices from fast5 output stored in S3")

    parser.add_argument('--database', '-d', action='store', dest='database', required=True, type=str,
                        help="location of simpleDB database file")
    parser.add_argument('--initialize', '-I', action='store_true', dest='initialize', required=False, default=False,
                        help="initialize database into file")
    parser.add_argument('--input_glob', '-i', action='store', dest='input', type=str, default=None,
                        help="glob for files to import into database")
    parser.add_argument('--stats', '-s', action='store_true', dest='stats', default=False,
                        help="print stats for the file")
    parser.add_argument('--read_ids', '-r', action='store', dest='read_ids', type=str, default=None,
                        help="file with one read_id per line. fast5 locations found in db will printed to stdout")


    args = parser.parse_args()
    return args

class Fast5LookupDB(object):
    def __init__(self, database_location, initialize=False):
        if not os.path.isfile(database_location) and not initialize:
            raise Exception("Database file {} does not exist. Perhaps it needs to be initialized"
                            .format(database_location))
        self.db = database_location
        if initialize:
            self._initialize_database()


    def _initialize_database(self):
        if os.path.isfile(self.db):
            print("Cannot initialize an existing database: {}".format(self.db), file=sys.stderr)
            return

        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            c.execute(
                '''CREATE TABLE fast5_root (
                    id INTEGER PRIMARY KEY, 
                    location TEXT UNIQUE
                )''' )
            c.execute(
                '''CREATE TABLE fast5_run (
                    id INTEGER PRIMARY KEY, 
                    run TEXT
                )''' )
            c.execute(
                '''CREATE TABLE read_id_to_fast5_location (
                    read_id TEXT PRIMARY KEY, 
                    fast5_location TEXT, 
                    fast5_root_id INTEGER,
                    run_id INTEGER,
                    FOREIGN KEY(fast5_root_id) REFERENCES fast5_root(id),
                    FOREIGN KEY(run_id) REFERENCES fast5_run(id)
                )''')
            conn.commit()

    def _get_or_insert_fast5_root(self, fast5_root):
        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            for row in c.execute('''SELECT id FROM fast5_root WHERE location=?''', (fast5_root,)):
                return row[0]
            result = c.execute('''INSERT INTO fast5_root(location) VALUES (?)''', (fast5_root,))
            conn.commit()
            return result.lastrowid

    def _get_or_insert_fast5_run(self, run_name):
        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            for row in c.execute('''SELECT id FROM fast5_run WHERE run=?''', (run_name,)):
                return row[0]
            result = c.execute('''INSERT INTO fast5_run(run) VALUES (?)''', (run_name,))
            conn.commit()
            return result.lastrowid

    @staticmethod
    def _get_read_query():
        return '''SELECT r.read_id, r.fast5_location, o.location, u.run 
                  FROM read_id_to_fast5_location r
                   LEFT JOIN fast5_root o ON r.fast5_root_id = o.id
                   LEFT JOIN fast5_run u on r.run_id = u.id
               '''

    @staticmethod
    def _db_read_conversion(db_read):
        return {
            READ_ID: db_read[0],
            FAST5_LOCATION: db_read[1] if db_read[2] is None else "{}{}".format(db_read[2], db_read[1]),
            RUN_NAME: db_read[3]
        }

    @staticmethod
    def _only_fast5_location(read):
        return read[FAST5_LOCATION]

    def insert_reads(self, reads, fast5_root=None):
        master_fast5_root_id = None if fast5_root is None else self._get_or_insert_fast5_root(fast5_root)
        run_map = { read[RUN_NAME] : self._get_or_insert_fast5_run(read[RUN_NAME]) for read in reads}
        def read_to_db_format(read):
            read_id, fast5_location, run_id, fast5_root_id = \
                read[READ_ID], read[FAST5_LOCATION], run_map[read[RUN_NAME]], None
            if fast5_root is not None and fast5_location.startswith(fast5_root):
                fast5_location = fast5_location.replace(fast5_root, '')
                fast5_root_id = master_fast5_root_id
            return (read[READ_ID], fast5_location, fast5_root_id, run_id,)
        db_reads = list(map(read_to_db_format, reads))

        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            result = c.executemany('INSERT OR REPLACE INTO read_id_to_fast5_location VALUES (?,?,?,?)', db_reads)
            conn.commit()
            return result.rowcount

    def get_all(self, only_fast5_location=False):
        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            results = c.execute(Fast5LookupDB._get_read_query())
            reads = list(map(Fast5LookupDB._db_read_conversion, results))
            if only_fast5_location:
                return list(map(Fast5LookupDB._only_fast5_location, reads))
            return reads

    def get_by_read_ids(self, read_ids, only_fast5_location=False):
        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            results = c.execute(
                "{} WHERE read_id in ({})".format(Fast5LookupDB._get_read_query(), ",".join(["?" for _ in read_ids])),
                read_ids)
            reads = list(map(Fast5LookupDB._db_read_conversion, results))
            if only_fast5_location:
                return list(map(Fast5LookupDB._only_fast5_location, reads))
            return reads

    def get_by_read_id(self, read_id, only_fast5_location=False):
        return self.get_by_read_ids([read_id], only_fast5_location)

    def get_stats(self):
        filesize = os.stat(self.db).st_size / (1 << 20)
        with closing(sqlite3.connect(self.db)) as conn:
            c = conn.cursor()
            read_count = 0
            for row in c.execute("SELECT COUNT(*) FROM read_id_to_fast5_location"):
                read_count = row[0]

            run_count = 0
            for row in c.execute("SELECT COUNT(*) FROM fast5_run"):
                run_count = row[0]

            root_count = 0
            for row in c.execute("SELECT COUNT(*) FROM fast5_root"):
                root_count = row[0]

            return read_count, run_count, root_count, filesize


def load_index_into_database(db, index_glob):
    import glob
    files = glob.glob(index_glob)
    files.sort()
    print("\nLoading reads from {} files matching {}".format(len(files), index_glob), file=sys.stderr)
    total_read_count = 0

    for file_idx, index_file in zip(range(len(files)), files):
        fast5_root, fast5_idx, read_id_idx, run_id_idx = None, None, None, None
        reads = list()
        with open(index_file, 'r') as input:
            for line in input:
                if line.startswith("##"):
                    if line.startswith("##{}:".format(FAST5_ROOT)):
                        fast5_root = line.strip()[len("##{}:".format(FAST5_ROOT)):]
                        # fast5_root = line.lstrip("##{}:".format(FAST5_ROOT))
                elif line.startswith("#"):
                    line = line.lstrip('#').split()
                    if fast5_idx is None:
                        for column, idx in zip(line, range(len(line))):
                            if column == FAST5_LOCATION: fast5_idx = idx
                            if column == READ_ID: read_id_idx = idx
                            if column == RUN_NAME: run_id_idx = idx
                        if None in [fast5_idx, read_id_idx, run_id_idx]:
                            raise Exception("Invalid header line in {}: {}".format(index_file, line))
                else:
                    line = line.split()
                    reads.append({FAST5_LOCATION:line[fast5_idx], READ_ID:line[read_id_idx], RUN_NAME:line[run_id_idx]})

        stored_count = db.insert_reads(reads, fast5_root)
        print("{}: Loaded {} / {} into database from {}".format(file_idx, stored_count, len(reads), index_file),
              file=sys.stderr)
        total_read_count += stored_count

    print("Loaded total of {} reads successfully.".format(total_read_count), file=sys.stderr)


def get_fast5_locations(db, read_filename, print_to_stdout=True):

    read_ids = list()
    with open(read_filename, 'r') as input:
        for line in input:
            if not line.startswith("#"): read_ids.append(line.strip())
    print("\nRead {} reads for downloading from {}".format(len(read_ids), read_filename), file=sys.stderr)

    fast5_locations = db.get_by_read_ids(read_ids, only_fast5_location=True)
    print("Found {} fast5 locations (of {} reads)".format(len(fast5_locations), len(read_ids)), file=sys.stderr)

    if print_to_stdout:
        for fast5 in fast5_locations:
            print(fast5)

    return fast5_locations

def main():
    args = parse_args()

    fast5_db = Fast5LookupDB(args.database, initialize=args.initialize)

    if args.input is not None:
        load_index_into_database(fast5_db, args.input)

    if args.stats:
        read_count, run_count, root_count, filesize = fast5_db.get_stats()
        print("{}:\n\tfile_size:  {} Mb\n\tread_count: {}\n\trun_count:  {}\n\troot_count: {}\n"
              .format(args.database, filesize, read_count, run_count, root_count))

    if args.read_ids is not None:
        get_fast5_locations(fast5_db, args.read_ids)




if __name__ == "__main__":
    main()