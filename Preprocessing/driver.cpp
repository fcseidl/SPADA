//
//  driver.cpp
//  Preprocessing
//
//  Created by Frank Seidl on 7/21/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#include <cassert>

#include "xcode_redirect.hpp"
#include "csvtools.hpp"
#include "topics.hpp"


using namespace std;

const int CORNER_SQUARE_SIZE = 4;

const string USAGE =
"Allowable usages:\n\
./preprocess.exe -wordbag [infile] [outfile] [bag_size]\n\
./preprocess.exe -corners [filename]\n\
./preprocess.exe -corners [filename] -ignore [ignore_rows]\n\
./preprocess.exe -ssf [small] [large] [small_out] [large_out]\n\
./preprocess.exe -join [small] [large] [out]\n\
./preprocess.exe -help\n";

enum class Mode { CORNERS, SSF, JOIN, BAG };

struct Arguments {
    Mode mode;
    int ignore_rows;         // only for corners mode
    int bag_size;            // only for wordbag mode
    vector<string> filenames;
}; // Arguments


void fail_if(bool condition) {
    if (condition) {
        cout << "Preprocessor is confused. The arguments provided must be in one of the forms below.\n";
        cout << USAGE;
        exit(1);
    } // if
} // fail_with_message()


Arguments process_args(int argc, char *argv[]) {
    Arguments result;
    for (int arg = 1; arg < argc; ++arg) {
        switch (argv[arg][1]) {
            case 'c':
                // corners mode, next arg must be a file name
                result.mode = Mode::CORNERS;
                fail_if(++arg >= argc);
                result.filenames.emplace_back(argv[arg]);
                break;
                
            case 'h':
                // help flag passed
                cout << USAGE;
                exit(0);
                
            case 'i':
                // ignore rows only in corners mode
                fail_if(result.mode != Mode::CORNERS);
                fail_if(++arg >= argc);
                result.ignore_rows = atoi(argv[arg]);
                break;
                
            case 'j':
                // join mode, next 3 args are file names
                result.mode = Mode::JOIN;
                while (++arg < argc) result.filenames.emplace_back(argv[arg]);
                fail_if(result.filenames.size() != 3);
                break;
                
            case 's':
                // SSF mode, next 4 args are file names
                result.mode = Mode::SSF;
                while (++arg < argc) result.filenames.emplace_back(argv[arg]);
                fail_if(result.filenames.size() != 4);
                break;
                
            case 'w':
                // word bag mode, next args are infile, outfile, bag size
                result.mode = Mode::BAG;
                while (++arg < argc - 1) result.filenames.emplace_back(argv[arg]);
                fail_if(result.filenames.size() != 2);
                result.bag_size = atoi(argv[arg]);
                break;
                
            default:
                // incorrect flag
                fail_if(true);
                break;
                
        } // switch
    } // for arg
    
    return result;
} // process_args()


// TODO: don't allow other delimiters besides ','
// TODO: debug corners mode

int main(int argc, char *argv[]) {
    // speeds up I/O
    ios_base::sync_with_stdio(false);

    // redirect input for Xcode debugging
    xcode_redirect(argc, argv);

    // get arguments
    Arguments args = process_args(argc, argv);

    switch (args.mode) {
        case Mode::CORNERS:
            print_corners(CORNER_SQUARE_SIZE, args.filenames[0], DEFAULT_DELIM, args.ignore_rows);
            break;

        case Mode::SSF: {
            csvJoiner joiner(args.filenames[0], args.filenames[1], DEFAULT_DELIM);
            joiner.sorted_shared_features(args.filenames[2], args.filenames[3]);
            break;
        }

        case Mode::JOIN: {
            csvJoiner joiner(args.filenames[0], args.filenames[1], DEFAULT_DELIM);
            joiner.join(args.filenames[2]);
            break;
        }
        
        case Mode::BAG: {
            WordBagger bagger(args.filenames[0], args.bag_size);
            bagger.write_csv(args.filenames[1], DEFAULT_DELIM);
            break;
        }

        default:
            assert(false);
            break;
    } // switch
    
    return 0;
} // main()
