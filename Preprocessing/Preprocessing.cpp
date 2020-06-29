//
//  Preprocessing.cpp
//  Preprocessing
//
//  Created by Frank Seidl on 6/24/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#include "Preprocessing.hpp"
#include "xcode_redirect.hpp"

using namespace std;


/*
 Return string containing first token of a line.
 */
std::string first_token(std::string line, const char * delims) {
    char * str = &(line[0]);
    char * token = strtok(str, delims);
    return string(token);
} // first_token()


/*
Take two csv files small_in and large_in. For fastest performance,
the file with fewer rows should be first. Write two new csv files,
small_out and large_out, containing rows of small_in and large_in
whose names (first tokens) appear in both files, in lexicographical
sorted order. Row names must be unique.

Runtime: linear in lengths of files, log-linear in number of shared
row names.
*/
void sorted_shared_features(const string &small_in,
                            const string &large_in,
                            const string &small_out,
                            const string &large_out,
                            const char * delims) {
    string row;
    unordered_map<string, string> small_rows_by_name;
    map<string, string> large_matched_rows;
    ifstream read_small(small_in), read_large(large_in);
    ofstream write_small(small_out), write_large(large_out);
    
    // read rows of first file, hash them by name
    while (getline(read_small, row)) {
        string row_name = first_token(row, delims);
        small_rows_by_name[row_name] = row;
    } // while reading 1
    
    // read rows of second file, storing them in map if they match rows in first
    while (getline(read_large, row)) {
        string row_name = first_token(row, delims);
        if (small_rows_by_name.count(row_name)) {
            large_matched_rows[row_name] = row;
        } // if
    } // while reading 2
    
    // write matching rows in sorted order into both outfiles
    for (auto &pair:large_matched_rows) {
        const string &row_name = pair.first;
        const string &large_row = pair.second;
        write_large << large_row << '\n';
        write_small << small_rows_by_name[row_name] << '\n';
    } // for matching rows
} // only_shared_features()


/*
Read a 2d array from a file. Rows are separated by newlines;
columns are separated by any delimiters which appear in the
parameter delims. Print the corner values of the array, using
the following format:

1 1 1 ... 3 1 3
4 1 3 ... 0 2 1
3 6 2 ... 0 2 0
...
6 3 0 ... 0 0 4
3 3 4 ... 1 3 2
1 1 3 ... 3 0 7

The size of the corner squares printed is determined by the parameter
num. If num is larger than half the number of rows or columns, then
entire rows or columns are printed without ellipses.

Finally, print the dimensions of the array.

Single-pass, O(num^2) space.
*/
void print_corners(int num, const string &filename, const char * delims) {
    string line;
    ifstream reader(filename);
    deque< deque<string> > left(1);
    deque< deque<string> > right(1);
    char * ptr;
    size_t row, col, num_cols, begin_right;
    
    // read first line, determine column index at which right corners begin
    getline(reader, line);
    num_cols = 0;
    ptr = strtok(&(line[0]), delims);
    while (ptr) {
        string token(ptr);
        if (num_cols < num) {   // left gets first num entries
            left[0].push_back(token);
        } else {                // rightmost num entries in right
            if (right[0].size() >= (size_t)num) right[0].pop_front();
            right[0].push_back(token);
        } // if-else
        ptr = strtok(nullptr, delims);
        ++num_cols;
    } // while token
    begin_right = num_cols - right[0].size();
    
    // iterate over successive rows starting at second row
    for (row = 1; getline(reader, line); ++row) {
        if (row == (size_t)num) {
            // if num rows read, then print top corners
            for (size_t r = 0; r < num; ++r) {
                print_tsv(left[r]);
                if (2 * num < num_cols) cout << "...\t";
                print_tsv(right[r]);
                cout << '\n';
            } // for r
            left.clear();
            right.clear();
        } else if (left.size() >= (size_t)num) {
            // if more than num rows remembered, discard earliest remembered row
            left.pop_front();
            right.pop_front();
        } // if-else
        
        left.push_back(deque<string>());
        right.push_back(deque<string>());
        
        // read new row
        col = 0;
        ptr = strtok(&(line[0]), delims);
        while (ptr) {
            string token(ptr);
            if (col < num) {   // left gets first num entries
                left.back().push_back(token);
            } else if (col >= begin_right) {
                right.back().push_back(token);
            } // if-else
            ptr = strtok(nullptr, delims);
            ++col;
        } // while token
    } // for row
    
    // print bottom corners
    if (row > 2 * num) {    // are we skipping middle rows?
        cout << "...\n";
    }
    for (size_t r = 0; r < left.size(); ++r) {
        print_tsv(left[r]);
        if (2 * num < num_cols) cout << "...\t";
        print_tsv(right[r]);
        cout << '\n';
    } // for r
    
    // print dimensions
    cout << "\nArray is " << row << " rows by " << num_cols << " columns.\n";
} // print_corners()


///*
// Read a 2d array from stdin, and print its first column to stdout.
// Columns delimiting characters are in the parameter delims.
// 
// Single-pass.
//*/
//void print_first_col(const char * delims) {
//    string line;
//    char * firsttok;
//    while (getline(cin, line)) {
//        firsttok = strtok(&(line[0]), delims);
//        cout << firsttok << '\n';
//    } // while
//} // print_first_col()


int main(int argc, char *argv[]) {
    // speeds up I/O
    ios_base::sync_with_stdio(false);
    
    // redirect input for Xcode debugging
    xcode_redirect(argc, argv);
    
    switch (argv[1][0]) {
        case 'c':   // corners
            print_corners(4, argv[2], DEFAULT_DELIMS);
            break;
        case 's':   // ssf
            sorted_shared_features(argv[2], argv[3], argv[4], argv[5], DEFAULT_DELIMS);
            break;
        default:
            cout << "Usage:\n"
                << "./preprocess.exe corners [filename]\n"
                << "or\n"
                << "./preprocess.exe ssf [small_in] [large_in] [small_out] [large_out]\n";
            return 1;
    } // switch
    return 0;
} // main
