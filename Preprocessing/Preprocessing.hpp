//
//  Preprocessing.hpp
//  Preprocessing
//
//  Created by Frank Seidl on 6/24/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#ifndef Preprocessing_hpp
#define Preprocessing_hpp

#include <stdio.h>
#include <iostream>     // cin, cout
#include <string>       // string
#include <deque>        // deque
#include <string.h>     // strtok
#include <algorithm>    // max
#include <iomanip>      // decimal precision

/*
 Print elements of a range separated by tabs.
 */
template<template<class...> class rangeType, class eltType>
void print_tsv(rangeType<eltType> & range) {
    for (eltType elt:range)
        std::cout << elt << '\t';
} // print_tsv


/*
 EFFECTS:
 Read a 2d array from stdin. Rows are separated by newlines;
 columns are separated by any delimiter. Print the array, omitting
 most values, using the following format:
 
 1 1 1 ... 3 1 3
 4 1 3 ... 0 2 1
 3 6 2 ... 0 2 0
 ...
 6 3 0 ... 0 0 4
 3 3 4 ... 1 3 2
 1 1 3 ... 3 0 7
 
 The size of the corner squares printed is determined by the parameter
 num.
 */
void print_corners(int num, char delim);

#endif /* Preprocessing_hpp */
