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
#include "removal.hpp"


using namespace std;

const char DEFAULT_DELIM = ',';

const size_t CORNER_SQUARE_SIZE = 4;

const string USAGE =
"Usage: please use one of the modes below.\n\
\n\
corners mode\n\
------------\n\
Supply the -corners flag, followed by the name of a rectangular csv file.\n\
The optional argument -ignore [k] will ignore the first k rows.\n\
\n\
sorted shared features mode\n\
---------------------------\n\
Supply the -ssf flag, followed by four filenames:\n\
-ssf [small] [large] [small_out] [large_out]\n\
\n\
join mode\n\
---------\n\
Supply the -join flag, followed by three filenames:\n\
-join [small] [large] [out]\n\
\n\
word bagging mode\n\
-----------------\n\
Supply the -wordbag flag, followed by two filenames and a whole number:\n\
-wordbag [infile] [outfile] [bag_size]\n\
\n\
removal mode\n\
------------\n\
Supply the -remove flag, followed by two file names and a string:\n\
-remove [infile] [outfile] [remove]\n\
\n\
alter mode\n\
----------\n\
Supply the -alter flag, followed by two file names and two strings:\n\
-alter [infile] [outfile] [remove] [replacements]\n\
new titles mode\n\
\n\
----------\n\
Supply the -newtitles flag, followed by three filenames and a whole number:\n\
-newtitles [infile] [outfile] [substfile] [col]\n\
\n\
help mode\n\
---------\n\
Supply the -help flag.\n\
\n\
For each mode, the delimiting character can be changed from a comma to\n\
a tab with the optional -tsv flag, or to a space with -psv. This flag must\n\
be passed first.\n";

enum class Mode { CORNERS, SSF, JOIN, BAG, REMOVE, ALTER, NEWTITLES };

struct Arguments {
    Mode mode;
    int ignore_rows;         // only for corners mode
    int bag_size;            // only for wordbag mode
    int col;                 // only for newtitles mode
    char delim;
    vector<string> strings;
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
    result.delim = DEFAULT_DELIM;
    for (int arg = 1; arg < argc; ++arg) {
        switch (argv[arg][1]) {
            case 'a':
                // alter mode, next 4 args are strings
                result.mode = Mode::ALTER;
                while (++arg < argc) result.strings.emplace_back(argv[arg]);
                fail_if(result.strings.size() != 4);
                break;
                
            case 'c':
                // corners mode, next arg must be a file name
                result.mode = Mode::CORNERS;
                fail_if(++arg >= argc);
                result.strings.emplace_back(argv[arg]);
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
                while (++arg < argc) result.strings.emplace_back(argv[arg]);
                fail_if(result.strings.size() != 3);
                break;
                
            case 'n':
                // newtitles mode, 3 filenames and one int
                result.mode = Mode::NEWTITLES;
                while (++arg < argc) result.strings.emplace_back(argv[arg]);
                fail_if(result.strings.size() != 4);
                result.col = atoi(result.strings.back().c_str());   // TODO: hacky
                break;
                
            case 'p':
                // -psv flag for space delimiter
                result.delim = ' ';
                break;
                
            case 'r':
                // REMOVE mode, next 3 args are strings
                result.mode = Mode::REMOVE;
                while (++arg < argc) result.strings.emplace_back(argv[arg]);
                fail_if(result.strings.size() != 3);
                break;
                
            case 's':
                // SSF mode, next 4 args are file names
                result.mode = Mode::SSF;
                while (++arg < argc) result.strings.emplace_back(argv[arg]);
                fail_if(result.strings.size() != 4);
                break;
            
            case 't':
                // -tsv flag for table delimiter
                result.delim = '\t';
                break;
                
            case 'w':
                // word bag mode, next args are infile, outfile, bag size
                result.mode = Mode::BAG;
                while (++arg < argc - 1) result.strings.emplace_back(argv[arg]);
                fail_if(result.strings.size() != 2);
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


// TODO: fully support other delimiters besides comma
// TODO: allow null tokens!

int main(int argc, char *argv[]) {
    // speeds up I/O
    ios_base::sync_with_stdio(false);

    // redirect input for Xcode debugging
    xcode_redirect(argc, argv);
    
    // get arguments
    Arguments args = process_args(argc, argv);

    switch (args.mode) {
        case Mode::ALTER:
            subst_chars(args.strings[0], args.strings[1], args.strings[2], args.strings[3]);
            break;
        
        case Mode::CORNERS:
            print_corners(CORNER_SQUARE_SIZE, args.strings[0], args.delim, args.ignore_rows);
            break;

        case Mode::SSF: {
            csvJoiner joiner(args.strings[0], args.strings[1], args.delim);
            joiner.sorted_shared_features(args.strings[2], args.strings[3]);
            break;
        }

        case Mode::JOIN: {
            csvJoiner joiner(args.strings[0], args.strings[1], args.delim);
            joiner.join(args.strings[2]);
            break;
        }
            
        case Mode::NEWTITLES:
            subst_title_col(args.strings[0], args.strings[1], args.strings[2], args.col, args.delim);
            break;

        case Mode::BAG: {
            WordBagger bagger(args.strings[0], args.bag_size);
            bagger.write_csv(args.strings[1], args.delim);
            break;
        }
        
        case Mode::REMOVE:
            remove_chars(args.strings[0], args.strings[1], args.strings[2]);
            break;

        default:
            assert(false);
            break;
    } // switch
    
    return 0;
} // main()
