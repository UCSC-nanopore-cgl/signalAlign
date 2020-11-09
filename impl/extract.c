//
// Created by Andrew Bailey on 10/25/18.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <eventAligner.h>

#define RAW_ROOT "/Raw/Reads/"

void usage() {
    fprintf(stderr, "\n\textract - pull fastqs from fast5 files \n\n");
    fprintf(stderr, "--help/-h: Display this super useful message and exit\n");
    fprintf(stderr, "-d: fast5 dir\n");
    fprintf(stderr, "-o: name of output fastq file\n");
    fprintf(stderr, "-r: search all immediate subdirectories\n");

}



int main(int argc, char *argv[]) {
    char *fast5dir = NULL;
    char *output_file = NULL;
    bool recursive = false;
    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"recursive",               no_argument,        0,  'r'},
                {"fast5dir",                required_argument,  0,  'd'},
                {"output",                  required_argument,  0,  'o'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:rd:o:",
                          long_options, &option_index);

        if (key == -1) {
//            usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 'd':
                fast5dir = stString_copy(optarg);
                break;
            case 'o':
                output_file = stString_copy(optarg);
                break;
            case 'r':
                recursive = true;
                break;
            default:
                usage();
                return 1;
        }
    }

    char* fastq_out_file= "";
    char* readdb_out_file = "";

    // sanity checks
    if (fast5dir == NULL) st_errAbort("--fast5dir/-d parameter is required\n");
    if (output_file == NULL) st_errAbort("--output_file/-o parameter is required\n");

    stList *find_dots = stString_splitByString(output_file, ".");

    if (stList_length(find_dots) == 1){
        fastq_out_file = stString_concat(output_file, ".fastq");
        readdb_out_file = stString_concat(output_file, ".readdb");

    } else {
        if (stString_eq(stList_get(find_dots, 1), "fastq") || stString_eq(stList_get(find_dots, 1), "fq")){
            readdb_out_file = stString_concat(stList_get(find_dots, 0), ".readdb");
            fastq_out_file = output_file;
        } else {
            st_errAbort("output file must have fastq or fq extension: %s\n", fastq_out_file);
        }
    }


    if (stFile_exists(fastq_out_file)){
        st_errAbort("Fastq file already exists: %s\n", fastq_out_file);
    }
    if (stFile_exists(readdb_out_file)){
        st_errAbort("readdb file already exists: %s\n", readdb_out_file);
    }
    if (recursive){
        char* fast5_subdir;
        char* fast5_subdir_path;
        stList* fast5_dirs = stFile_getFileNamesInDirectory(fast5dir);
        for (int i = 0; i < stList_length(fast5_dirs); i++){
            fast5_subdir = stList_get(fast5_dirs, i);
            fast5_subdir_path = path_join_two_strings(fast5dir, fast5_subdir);
            if (stFile_isDir(fast5_subdir_path)){
                printf("%s\n", fast5_subdir_path);
                write_fastq_and_readdb_file1(fast5_subdir_path, fastq_out_file, readdb_out_file);
            }
        }
    } else {
        write_fastq_and_readdb_file1(fast5dir, fastq_out_file, readdb_out_file);
    }
    free(fast5dir);
    free(output_file);
    free(fastq_out_file);
    free(readdb_out_file);
    stList_destruct(find_dots);
}
