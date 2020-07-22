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
#include <vector>
#include <string.h>         // strtok
#include <algorithm>        // max
#include <iomanip>          // decimal precision
#include <utility>          // pair
#include <fstream>          // ifstream, ofstream
#include <map>
#include <unordered_map>


const char DEFAULT_DELIM = ',';


/*
 Return string containing first token of a line.
 */
std::string first_token(const std::string &line, char delims);


// Combines data from two different csv files.
class csvJoiner {
public:
    
    /*
     Create a csvJoiner of two csv files.
     */
    csvJoiner(std::string &small,
              std::string &large,
              char delim);
    
    /*
     Write a new csv file, out, with a row for each row name that is
     in both small and large. Each row is a concatenation of the corresponding
     rows in small and large. The first row of the file contains two values,
     [number of samples in small], [number of samples in large].
     */
    void join(std::string &out) const;
    
    /*
     Write two new csv files, small_out and large_out, containing rows of
     small and large whose names (first tokens) appear in both files,
     in lexicographical sorted order. Row names must be unique.
     
     Runtime: linear in lengths of files, log-linear in number of shared
     row names.
     */
    void sorted_shared_features(const std::string &small_out,
                                const std::string &large_out) const;
    
private:
    
    char delim;
    
    // header rows
    std::string small_header, large_header;
    
    // map row names to rows
    std::unordered_map<std::string, std::string> rows_hash;
    std::map<std::string, std::string> rows_bst;
}; // csvJoiner


/*
 Print elements of a range separated by tabs.
 */
template<template<class...> class rangeType, class eltType>
void print_tsv(const rangeType<eltType> & range) {
    for (eltType elt:range)
        std::cout << elt << '\t';
} // print_tsv


/*
 Read a 2d array from a file. Rows are separated by newlines;
 columns are separated by a delimiting character. Print the corner
 values of the array, using the following format:

 1 1 1 ... 3 1 3
 4 1 3 ... 0 2 1
 3 6 2 ... 0 2 0
 ...
 6 3 0 ... 0 0 4
 3 3 4 ... 1 3 2
 1 1 3 ... 3 0 7
 
 The size of the corner squares printed is determined by the parameter
 num. If num is larger than half the number of rows or columns, then
 entire rows or columns are printed without ellipses. Optionally, the
 first rows can be ignored.
 
 Finally, print the dimensions of the array.
 
 Single-pass, O(num^2) space.
*/
void print_corners(int num,
                   const std::string &filename,
                   char delim,
                   int ignore_rows = 0);

#endif /* Preprocessing_hpp */
