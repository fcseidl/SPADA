//
//  removal.cpp
//  Preprocessing
//
//  Created by Frank Seidl on 7/27/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#include "removal.hpp"

using namespace std;

/*
 Write contents of one file to a new file, excluding every character which appears in a string.
 */
void remove_chars(string& infile, string& outfile, string& chars) {
    string line;
    ifstream reader(infile);
    ofstream writer(outfile);
    while (getline(reader, line)) {
        for (char c:line) {
            if (chars.find(c) == string::npos) writer << c;
        }
        writer << '\n';
    } // while
} // remove_chars()
