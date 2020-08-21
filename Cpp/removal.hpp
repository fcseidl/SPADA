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
#include <cassert>

/*
 Write contents of one file to a new file, excluding every character which appears in a string.
 Second file must be different!
 */
void remove_chars(const std::string& infile,
                  const std::string& outfile,
                  const std::string& chars);


/*
 Write contents of one file to a new file, performing character substitutions according to
 in[i] -> out[i]. Second file must be different!
 */
void subst_chars(const std::string& infile,
                 const std::string& outfile,
                 const std::string& in,
                 const std::string& out);

#endif /* removal_hpp */
