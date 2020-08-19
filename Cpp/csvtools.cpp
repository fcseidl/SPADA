//
//  csvtools.cpp
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
    string row, row_name;
    ifstream read_small(small), read_large(large);
    
    // read header rows
    getline(read_small, small_header);
    getline(read_large, large_header);
    
    // read rows of first file, hash them by name
    while (getline(read_small, row)) {
        row_name = first_token(row, delim);
        rows_hash[row_name] = row;
    } // while reading 1
    
    // read rows of second file, storing them in map if they match rows in first
    while (getline(read_large, row)) {
        row_name = first_token(row, delim);
        if (rows_hash.count(row_name)) {
            rows_bst[row_name] = row;
        } // if
    } // while reading 2
} // ctor

/*
 Write a new csv file, out, with a row for each row name that is
 in both small and large. Each row is a concatenation of the corresponding
 rows in small and large. The first row of the file contains two values,
 [number of samples in small], [number of samples in large].
 
 // TODO: support non comma delims?
 // TODO: avoid inserting consecutive delimiting characters
 */
void csvJoiner::join(std::string &out) const {
    ofstream write_out(out);
    
    // determine number of samples--each is preceded by a comma
    const string &small_row = rows_hash.begin()->second;
    const string &large_row = rows_bst.begin()->second;
    size_t small_width = (size_t)count(small_row.begin(), small_row.end(), delim);
    size_t large_width = (size_t)count(large_row.begin(), large_row.end(), delim);
    write_out << small_width << delim << large_width << '\n';
    
    // write header row
    write_out << small_header;
    string trimmed = large_header.substr(large_header.find(delim));
    write_out << trimmed << '\n';
    
    // write concatenated rows
    for (auto &pair:rows_bst) {
        const string &row_name = pair.first;
        const string &large_row = pair.second;
        const string &small_row = rows_hash.at(row_name);
        // must trim out name column of large file
        trimmed = large_row.substr(large_row.find(delim));
        write_out << small_row << delim << trimmed << '\n';
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
    
    // write headers
    write_small << small_header << '\n';
    write_large << large_header << '\n';
    
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
columns are separated by separated by a delimiting character.
Print the corner values of the array, using the following format:

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
void print_corners(size_t num,
                   const string &filename,
                   char delim,
                   int ignore_rows) {
    string line;
    char * token;
    ifstream readfile(filename);
    deque< deque<string> > left;
    deque< deque<string> > right;
    size_t row = 0;             // row index
    size_t width;               // total number of columns
    size_t begin_right;         // column index of beginning of right-side corners
    
    // ignore first rows
    for (int i = 0; i < ignore_rows; ++i) getline(readfile, line);
    
    // read first line, determine width and begin_right
    getline(readfile, line);
    width = (size_t)count(line.begin(), line.end(), delim) + 1;
    begin_right = max(width - num, num);
    
    // iterate over successive rows
    do {
        if ((size_t)count(line.begin(), line.end(), delim) + 1 != width) {
            cout << "\nError: cannot print corners of ragged array.\n";
            exit(1);
        }
        // read ends of current row
        string delim_str = "_";
        delim_str[0] = delim;
        token = strtok(&(line[0]), delim_str.c_str());
        left.push_back(deque<string>());
        right.push_back(deque<string>());
        for (size_t col = 0; col < width; ++col) {
            if (col < num) {
                left.back().emplace_back(token);
            } else if (col >= begin_right) {
                right.back().emplace_back(token);
            } // if-else
            token = strtok(nullptr, delim_str.c_str());
        } // for col
        
        // print if num rows reached, discard oldest row if past num rows
        if (row == num) {
            for (size_t r = 0; r < num; ++r) {
                print_tsv(left[r]);
                if (2 * num < width) cout << "...\t";
                print_tsv(right[r]);
                cout << '\n';
            } // for r
            left.clear();
            right.clear();
        } else if (left.size() > num) {
            left.pop_front();
            right.pop_front();
        }
        
        // proceed to next row
        ++row;
    } while (getline(readfile, line));
    
    // print "..." if rows were skipped, then bottom corners
    if (2 * num < row) cout << "...\n";
    for (size_t r = 0; r < left.size(); ++r) {
        print_tsv(left[r]);
        if (2 * num < width) cout << "...\t";
        print_tsv(right[r]);
        cout << '\n';
    } // for r
    
    // print dimensions
    cout << "\nArray is " << row << " rows by " << width << " columns.\n";
} // print_corners()

