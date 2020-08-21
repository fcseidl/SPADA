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
 Second file must be different!
 */
void remove_chars(const string& infile,
                  const string& outfile,
                  const string& chars) {
    string line;
    ifstream reader(infile);
    ofstream writer(outfile);
    while (getline(reader, line)) {
        for (char c:line) {
            if (chars.find(c) == string::npos) writer << c;
        } // for c
        writer << '\n';
    } // while
} // remove_chars()


/*
 Write contents of one file to a new file, performing character substitutions according to
 in[i] -> out[i]. Second file must be different!
 */
void subst_chars(const string& infile,
                 const string& outfile,
                 const string& in,
                 const string& out) {
    assert(in.size() == out.size());
    size_t pos;
    string line;
    ifstream reader(infile);
    ofstream writer(outfile);
    while (getline(reader, line)) {
        for (char c:line) {
            if ((pos = in.find(c)) == string::npos) {
                writer << c;
            } else {
                writer << out[pos];
            } // if-else
        } // for c
        writer << '\n';
    } // while
} // subst_chars()
