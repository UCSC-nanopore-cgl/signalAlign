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
    fprintf(stderr, "-o: output file\n");
}



int main(int argc, char *argv[]) {
    char *fast5dir = NULL;
    char *out_file = NULL;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"fast5dir",                required_argument,  0,  'd'},
                {"out_file",                required_argument,  0,  'o'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:d:o:",
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
                out_file = stString_copy(optarg);
                break;
            default:
                usage();
                return 1;
        }
    }

    // sanity checks
    if (fast5dir == NULL) st_errAbort("--fast5dir/-d parameter is required\n");
    if (out_file == NULL) st_errAbort("--out_file/-o parameter is required\n");
    if (stFile_exists(out_file)){
        st_errAbort("--out_file/-o parameter cannot already exist\n");
    }
    write_fastqs_to_file(fast5dir, out_file);

}
