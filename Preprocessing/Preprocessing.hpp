//
//  Preprocessing.hpp
//  Preprocessing
//
//  Created by Frank Seidl on 6/24/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#ifndef Preprocessing_hpp
#define Preprocessing_hpp

#include <iostream>         // cin, cout
#include <string>
#include <deque>
#include <string.h>         // strtok
#include <algorithm>        // max
#include <iomanip>          // decimal precision
#include <utility>          // pair
#include <fstream>          // ifstream, ofstream
#include <map>
#include <unordered_map>


/*
 Return string containing first token of a line.
 */
std::string first_token(std::string line, const char * delims);


/*
 Take two csv files small_in and large_in. For fastest performance,
 the file with fewer rows should be first. Write two new csv files,
 small_out and large_out, containing rows of small_in and large_in
 whose names (first tokens) appear in both files, in lexicographical
 sorted order. Row names must be unique.
 
 Runtime: linear in lengths of files, log-linear in number of shared
 row names.
 */
void sorted_shared_features(const std::string &small_in,
                            const std::string &small_out,
                            const std::string &large_in,
                            const std::string &large_out,
                            const char * delims);


/*
 Print elements of a range separated by tabs.
 */
template<template<class...> class rangeType, class eltType>
void print_tsv(const rangeType<eltType> & range) {
    for (eltType elt:range)
        std::cout << elt << '\t';
} // print_tsv


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

Finally, print the dimensions of the array.

Single-pass, O(num^2) space.
*/
void print_corners(int num, char delim);


/*
 Read a 2d array from stdin, and print its first column to stdout.
 Columns delimiting characters are in the parameter delims.
 
 Single-pass.
*/
void print_first_col();

#endif /* Preprocessing_hpp */
