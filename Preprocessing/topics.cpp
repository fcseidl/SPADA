//
//  topics.cpp
//  Preprocessing
//
//  Created by Frank Seidl on 7/21/20.
//  Copyright © 2020 Frank Seidl. All rights reserved.
//

#include "topics.hpp"

using namespace std;

/*
 Convert uppercase to lowercase, remove nonalphabetical characters.
 If this produces an empty string, replace with "_".
 */
void only_lowercase(string& str) {
    not_a_letter pred;
    remove_if(str.begin(), str.end(), pred);
    for (char& c:str) {
        if ('A' <= c && c <= 'Z') {
            c = c - 'A' + 'a';
        } else if (pred(c)) {
            c = '\0';   // terminate string after all letters
        }
    } // for c
    // empty string?
    if (str == "") str = "_";
} // only_lowercase()


//------------------------WordBagger implementation---------------------------//

/*
 Create a WordBagger for a file. Reading file occurs during construction.
 */
WordBagger::WordBagger(string& infile, int bag_size_in)
: bag_size(bag_size_in), num_bags(0), reader(infile) {
    // read bags until none reamin
    while (read_bag()) {}
} // ctor

/*
 Write word count data into a csv-like format with a specified delimiting
 character.
 */
void WordBagger::write_csv(string& outfile, char delim) const {
    ofstream write_out(outfile);
    
    // write top column, with number label of each bag
    write_out << "word/bag" << delim;
    for (size_t label = 0; label < num_bags; ++label)
        write_out << label << delim;
    write_out << '\n';
    
    // write bag counts for each word in no particular order
    for (auto& tuple:bag_counts) {
        // write the word as row header
        write_out << tuple.first << delim;
        
        // write count of this word in each bag
        const vector<pair<size_t, size_t>>& word_counts = tuple.second;
        size_t bag = 0;
        size_t index = 0;   // index into word_counts
        while (bag < num_bags) {
            if (word_counts[index].first == bag) {
                write_out << word_counts[index].second;
                ++index;
            } else {
                write_out << 0;
            } // if-else
            write_out << delim;
            ++bag;
        } // for each bag
        write_out << '\n';
    } // for each word
} // write_csv()

/*
 Read one bag, return false if EOF was reached.
 */
bool WordBagger::read_bag() {
    string word;
    for (int i = 0; i < bag_size; ++i) {
        if (reader >> word) {
            only_lowercase(word);
            vector<pair<size_t, size_t>>& word_counts = bag_counts[word];
            if (!word_counts.empty() && word_counts.back().first == num_bags) {
                // this word has already appeared in the current bag
                word_counts.back().second++;
            } else {
                // this is the word's first appearance in current bag
                word_counts.emplace_back(num_bags, 1);
            } // if-else word already appeared
        } else {
            // end of file reached
            return false;
        } // if-else file not fully read
    } // for words
    // still more to read
    ++num_bags;
    return true;
} // read_bag()

//----------------------end WordBagger implementation-------------------------//
