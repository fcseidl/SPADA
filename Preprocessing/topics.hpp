//
//  topics.hpp
//  Preprocessing
//
//  Created by Frank Seidl on 7/21/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#ifndef topics_hpp
#define topics_hpp

#include <string>
#include <fstream>          // ifstream, ofstream
#include <algorithm>        // remove_if
#include <unordered_map>
#include <vector>
#include <utility>          // pair


// predicate to identify non-alphabetical characters
struct not_a_letter {
    bool operator()(char c) {
        if ('a' <= c && c <= 'z') return false;
        if ('A' <= c && c <= 'Z') return false;
        return true;
    }
};

/*
Convert uppercase to lowercase, remove nonalphabetical characters.
If this produces an empty string, replace with "_".
*/
void only_lowercase(std::string& str);


/*
 Class which reads blocks of words from a text file, and writes their
 "bag of words" vector representations to a csv file.
 
 STEPS:
 1. store data in a hash table where words map to lists of (block number, count) pairs
 2. print the data one row at a time
*/
class WordBagger {
public:
    
    /*
     Create a WordBagger for a file. Reading file occurs during construction.
     */
    WordBagger(std::string &infile, int bag_size_in);
    
    /*
     Write word count data into a csv-like format with a specified delimiting
     character.
     */
    void write_csv(std::string &outfile, char delim) const;
    
private:
    
    const int bag_size;
    size_t num_bags;
    std::ifstream reader;
    std::unordered_map< std::string, std::vector<std::pair<size_t, size_t>> > bag_counts;
    
    /*
     Read one bag, return false if EOF was reached.
    */
    bool read_bag();
    
}; // WordBagger


#endif /* topics_hpp */
