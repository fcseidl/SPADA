//
//  Preprocessing.cpp
//  Preprocessing
//
//  Created by Frank Seidl on 6/24/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#include "csvtools.hpp"
#include "xcode_redirect.hpp"

using namespace std;


/*
 Return string containing first token of a line.
 */
string first_token(const string &line, char delim) {
    size_t len = line.find(delim);
    return line.substr(0, len);
} // first_token()


//------------------------csvJoiner implementation---------------------------//

/*
 Create a csvJoiner of two csv files.
 */
csvJoiner::csvJoiner(string &small,
                     string &large,
                     char delim) : delim(delim) {
    string row;
    ifstream read_small(small), read_large(large);
    
    // read rows of first file, hash them by name
    while (getline(read_small, row)) {
        string row_name = first_token(row, delim);
        rows_hash[row_name] = row;
    } // while reading 1
    
    // read rows of second file, storing them in map if they match rows in first
    while (getline(read_large, row)) {
        string row_name = first_token(row, delim);
        if (rows_hash.count(row_name)) {
            rows_bst[row_name] = row;
        } // if
    } // while reading 2
} // ctor

/*
 Write a new csv file, out, with a row for each row name that is
 in both small and large. Each row is a concatenation of the corresponding
 rows in small andlarge. The first row of the file contains two values,
 [number of samples in small], [number of samples in large].
 */
void csvJoiner::join(std::string &out) const {
    ofstream write_out(out);
    
    // determine number of samples--each is preceded by a comma
    const string &small_row = rows_hash.begin()->second;
    const string &large_row = rows_bst.begin()->second;
    size_t small_width = (size_t)count(small_row.begin(), small_row.end(), delim);
    size_t large_width = (size_t)count(large_row.begin(), large_row.end(), delim);
    write_out << small_width << ',' << large_width << '\n';
    
    // write concatenated rows
    for (auto &pair:rows_bst) {
        const string &row_name = pair.first;
        const string &large_row = pair.second;
        const string &small_row = rows_hash.at(row_name);
        // must trim out name column of large file
        string trimmed = large_row.substr(large_row.find(delim));
        write_out << small_row << ',' << trimmed << '\n';
    } // for matching rows
} // join()

/*
 Write two new csv files, small_out and large_out, containing rows of
 small and large whose names (first tokens) appear in both files,
 in lexicographical sorted order. Row names must be unique.
 
 Runtime: linear in lengths of files, log-linear in number of shared
 row names.
 */
void csvJoiner::sorted_shared_features(const string &small_out,
                                       const string &large_out) const {
    ofstream write_small(small_out), write_large(large_out);
    
    // write matching rows in sorted order into both outfiles
    for (auto &pair:rows_bst) {
        const string &row_name = pair.first;
        const string &large_row = pair.second;
        write_large << large_row << '\n';
        write_small << rows_hash.at(row_name) << '\n';
    } // for matching rows
} // sorted_shared_features()

//----------------------end csvJoiner implementation-------------------------//


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
// TODO: this is some spaghetti code right here, oh boy
void print_corners(int num,
                   const string &filename,
                   char delim,
                   int ignore_rows) {
    string line;
    ifstream reader(filename);
    deque< deque<string> > left(1);
    deque< deque<string> > right(1);
    char * ptr;
    size_t row, col, num_cols, begin_right;
    
    // ignore first rows
    for (int i = 0; i < ignore_rows; ++i) getline(reader, line);
    
    // read first line, determine column index at which right corners begin
    getline(reader, line);
    num_cols = 0;
    ptr = strtok(&(line[0]), &delim);
    while (ptr) {
        string token(ptr);
        if (num_cols < (size_t)num) {   // left gets first num entries
            left[0].push_back(token);
        } else {                // rightmost num entries in right
            if (right[0].size() >= (size_t)num) right[0].pop_front();
            right[0].push_back(token);
        } // if-else
        ptr = strtok(nullptr, &delim);
        ++num_cols;
    } // while token
    begin_right = num_cols - right[0].size();
    
    // iterate over successive rows starting at second row
    for (row = 1; getline(reader, line); ++row) {
        if (row == (size_t)num) {
            // if num rows read, then print top corners
            for (size_t r = 0; r < (size_t)num; ++r) {
                print_tsv(left[r]);
                if (2 * (size_t)num < num_cols) cout << "...\t";
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
        ptr = strtok(&(line[0]), &delim);
        while (ptr) {
            string token(ptr);
            if (col < (size_t)num) {   // left gets first num entries
                left.back().push_back(token);
            } else if (col >= begin_right) {
                right.back().push_back(token);
            } // if-else
            ptr = strtok(nullptr, &delim);
            ++col;
        } // while token
    } // for row
    
    // print bottom corners
    if (row > 2 * (size_t)num) {    // are we skipping middle rows?
        cout << "...\n";
    }
    for (size_t r = 0; r < left.size(); ++r) {
        print_tsv(left[r]);
        if (2 * (size_t)num < num_cols) cout << "...\t";
        print_tsv(right[r]);
        cout << '\n';
    } // for r
    
    // print dimensions
    cout << "\nArray is " << row << " rows by " << num_cols << " columns.\n";
} // print_corners()

