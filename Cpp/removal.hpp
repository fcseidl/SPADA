//
//  removal.hpp
//  Preprocessing
//
//  Created by Frank Seidl on 7/27/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#ifndef removal_hpp
#define removal_hpp

#include <string>
#include <fstream>       // ifstream, ofstream
#include <algorithm>     // remove_if

//// functor to determine whether a character belongs to a string
//struct in_string {
//    const std::string& str;
//    bool operator()(char c) { return str.find(c) < str.size(); }
//};

/*
 Write contents of one file to a new file, excluding every character which appears in a string.
 */
void remove_chars(std::string& infile, std::string& outfile, std::string& chars);

#endif /* removal_hpp */
