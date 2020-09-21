#include <getopt.h>
#include "signalMachineUtils.h"
#include "eventAligner.h"


void usage() {
    fprintf(stderr, "\n\tsignalMachine - Align ONT ionic current to a reference sequence\n\n");
    fprintf(stderr, "--help: Display this super useful message and exit\n");
    fprintf(stderr, "-f: fast5 file\n");
    fprintf(stderr, "-m: model file\n");
    fprintf(stderr, "-p: path to save events in Fast5\n");
    fprintf(stderr, "-N: raw nucleotide string\n");
    fprintf(stderr, "-n: file with nucleotide string\n");
    fprintf(stderr, "-r: set read as rna \n");
    fprintf(stderr, "    (exactly one of -N or -n must be set)\n");
    fprintf(stderr, "-w: write failed alignments\n");
}


int main(int argc, char *argv[]) {
    char *fast5File = NULL;
    char *modelFile = NULL;
    char *savePath = NULL;
    char *nucleotides = NULL;
    char *nucleotideFile = NULL;
    bool writeFailedAlignment = false;
    bool rna = false;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"fast5",                   required_argument,  0,  'f'},
                {"model",                   required_argument,  0,  'm'},
                {"save_path",               required_argument,  0,  'p'},
                {"nucleotides",             required_argument,  0,  'N'},
                {"nucleotide_file",         required_argument,  0,  'n'},
                {"force_write",             no_argument,        0,  'w'},
                {"rna",                     no_argument,        0,  'r'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:f:m:p:N:n:wr",
                          long_options, &option_index);

        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 'f':
                fast5File = stString_copy(optarg);
                break;
            case 'm':
                modelFile = stString_copy(optarg);
                break;
            case 'p':
                savePath = stString_copy(optarg);
                break;
            case 'N':
                nucleotides = stString_copy(optarg);
                break;
            case 'n':
                nucleotideFile = stString_copy(optarg);
                break;
            case 'r':
                rna = true;
                break;
            case 'w':
                writeFailedAlignment = true;
                break;
            default:
                usage();
                return 1;
        }
    }

    // sanity checks
    if (fast5File == NULL) st_errAbort("--fast5/-f parameter is required\n");
    if (modelFile == NULL) st_errAbort("--model/-m parameter is required\n");
    if (savePath == NULL)  st_errAbort("--save_path/-p parameter is required\n");
    if ((nucleotides == NULL ? 1 : 0) + (nucleotideFile == NULL ? 1 : 0) != 1) {
      st_errAbort("exactly one of --nucleotides/-N and --nucleotide_file/-n parameter is required\n");
    }
    // do we need to get this from the fasta file?  it better just be one entry
    if (nucleotides == NULL) {
        FILE* fastaFile = fopen(nucleotideFile, "r");

        struct List *seqs = constructZeroLengthList(1, NULL);
        struct List *seqLengths = constructZeroLengthList(1, NULL);
        struct List *fastaNames = constructZeroLengthList(1, NULL);

        fastaRead(fastaFile, seqs, seqLengths, fastaNames);

        nucleotides = (char*) listRemoveFirst(seqs);
    }

    // do the alignment and return the status code
    return load_from_raw2(fast5File, modelFile, nucleotides, savePath, writeFailedAlignment, rna);
}
