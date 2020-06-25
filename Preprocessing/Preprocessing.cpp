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
 Read a 2d array from stdin. Rows are separated by newlines;
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
 
 Single-pass, O(num^2) space.
 */
void print_corners(int num, const char * delims) {
    string line;
    deque< deque<string> > left(1);
    deque< deque<string> > right(1);
    char * ptr;
    size_t row, col, num_cols, begin_right;
    
    // read first line, determine column index at which right corners begin
    getline(cin, line);
    num_cols = 0;
    ptr = strtok(&(line[0]), delims);
    while (ptr) {
        string token(ptr);
        token = token.substr(0, 3);
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
    for (row = 1; getline(cin, line); ++row) {
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
            token = token.substr(0, 3);
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
} // print_corners()


int main(int argc, char *argv[]) {
    // speeds up I/O
    ios_base::sync_with_stdio(false);
    
    // redirect input for Xcode debugging
    xcode_redirect(argc, argv);
    
    // format decimal numbers
    cout << setprecision(2) << fixed;

    print_corners(5, " \t");
    return 0;
} // main
