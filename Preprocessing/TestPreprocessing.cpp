//
//  TestPreprocessing.cpp
//  Preprocessing
//
//  Created by Frank Seidl on 6/24/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//
//  Simple test suite.
//

#include <stdio.h>
#include "Preprocessing.hpp"

using namespace std;

void test_print_tsv() {
    cout << "...testing print_tsv()...\n"
        << "expect: 1\t2\t3\t4\t\n"
        << "receive: ";
    deque<int> nums = { 1, 2, 3, 4 };
    print_tsv(nums);
    cout << '\n';
    
    cout << "expect: green\teggs\tand\tthe\tabsurd\n"
    << "receive: ";
    deque<const char*> words = {
        "green",
        "eggs",
        "and",
        "the",
        "absurd"
    };
    print_tsv(words);
    cout << '\n';
}

int main() {
    test_print_tsv();
    return 0;
}


