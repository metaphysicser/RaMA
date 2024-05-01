/*
 * Copyright [2024] [MALABZ_UESTC Pinglu Zhang]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

 // Author: Pinglu Zhang
 // Contact: pingluzhang@outlook.com
 // Created: 2024-02-20
#pragma once

#include "logging.h"
#include "gsacak.h"

#include <deque>
#include <map>
#include <cmath> 

// Structure to represent a rare match, including counts of occurrences in the first and second sequences,
// the length of the match, positions of the match, and the type of positions (first or second sequence).
struct RareMatch {
    uint_t first_count; // Number of occurrences in the first sequence.
    uint_t second_count; // Number of occurrences in the second sequence.
    uint_t match_length; // Length of the match.
    std::vector<uint_t> match_pos; // Positions of the match.
    std::vector<bool> pos_type; // Indicates the sequence of each position (0 for first, 1 for second).

    uint_t min_key; // Minimum key for sorting and identifying unique matches.

    // Default constructor initializes counts and min_key to zero and max value, respectively.
    RareMatch() : first_count(0), second_count(0), match_length(0), min_key(U_MAX) {}

    // Parameterized constructor to initialize the match with specific details.
    RareMatch(uint_t _match_length, const std::vector<uint_t> _match_pos, const std::vector<bool> _pos_type)
        : match_length(_match_length), match_pos(_match_pos), pos_type(_pos_type) {
        // Calculate the minimum key based on position and match length.
        min_key = *std::min_element(match_pos.begin(), match_pos.end()) + match_length;
        // Count occurrences in each sequence.
        first_count = std::count(pos_type.begin(), pos_type.end(), false);
        second_count = pos_type.size() - first_count;
    }
};

// Alias for a vector of RareMatch objects.
using RareMatches = std::vector<RareMatch>;
// Map to associate a unique key with each RareMatch.
using RareMatchMap = std::map<uint_t, RareMatch>;

// Structure to represent a pair of matching positions between sequences,
// including their starting positions, match length, and a weight for scoring.
struct RareMatchPair {
    uint_t first_pos; // Starting position in the first sequence.
    uint_t second_pos; // Starting position in the second sequence.
    uint_t match_length; // Length of the match.
    double weight; // Weight for scoring the match pair.

    // Operator overloading for sorting based on first and second positions.
    bool operator<(const RareMatchPair& other) const {
        if (first_pos != other.first_pos) {
            return first_pos < other.first_pos;
        }
        return second_pos < other.second_pos;
    }

    // Helper function to check if two RareMatchPairs are seamlessly adjacent
    bool isAdjacent(const RareMatchPair& next) const {
        return (first_pos + match_length == next.first_pos) &&
            (second_pos + match_length == next.second_pos);
    }

    // Helper function to check overlap
    bool hasOverlap(const RareMatchPair& next) const {
        return (first_pos + match_length > next.first_pos) ||
            (second_pos + match_length > next.second_pos);
    }

    // Helper function to merge two RareMatchPairs
    void mergeWith(const RareMatchPair& next) {
        // Assuming the match_length is the only attribute that needs to be updated
        match_length += next.match_length;
    }
};

// Alias for a vector of RareMatchPair objects.
using RareMatchPairs = std::vector<RareMatchPair>;

void saveRareMatchPairsToCSV(const RareMatchPairs& pairs, const std::string& filename, uint_t fst_len);
RareMatchPairs readRareMatchPairsFromCSV(const std::string& filename, uint_t fst_len);

// Represents an interval within the LCP (Longest Common Prefix) array.
class LCPInterval {
private:
    const std::vector<int_t>& LCP; // Reference to the LCP array.
    uint_t interval_size; // Size of the interval.

    uint_t left; // Left boundary of the interval.
    uint_t right; // Right boundary of the interval.

    std::deque<uint_t> min_deque; // Deque to maintain the minimum LCP value within the interval.
    uint_t min_LCP_value; // Cached minimum LCP value in the current interval.

public:
    // Constructor initializes the interval with a reference to the LCP array and interval size.
    explicit LCPInterval(const std::vector<int_t>& LCP_array, uint_t interval_size);

    // Moves the interval one position to the right and updates the minimum LCP value.
    void slideRight();

    // Returns the minimum LCP value within the interval.
    int_t getMinLCP() const;

    // Checks if the current interval represents a rare interval based on specific criteria.
    bool isRareInterval() const;

    // Checks if the right boundary of the interval has reached the end of the LCP array.
    bool isRightAtEnd() const;

    // Returns the boundaries of the current interval.
    std::pair<uint_t, uint_t> getboundary() const;
};


// Finds rare matches within concatenated sequences using LCP array.
class RareMatchFinder {
private:
    unsigned char* concat_data; // Concatenated sequence data.

    std::vector<uint_t> SA; // Suffix Array.
    std::vector<int_t> LCP; // Longest Common Prefix array.
    std::vector<int_da> DA; // Document Array indicating sequence origin.

    uint_t first_seq_start;
    uint_t first_seq_len; // Length of the first sequence.

    uint_t second_seq_start;
    uint_t second_seq_len; // Length of the second sequence.

    uint_t min_seq_len; // Minimum length of the two sequences.
    uint_t concat_seq_len; // Total length of the concatenated sequence.

    // Retrieves match positions and types from a given boundary.
    void getMatchPosAndType(std::pair<uint_t, uint_t> boundary, std::vector<uint_t>& match_pos, std::vector<bool>& pos_type);

    // Expands rare matches to the left within the rare match map.
    void leftExpandRareMatchMap(RareMatchMap& rare_match_map);

    // Expands match positions to the left and returns the number of expanded positions.
    uint_t leftExpand(std::vector<uint_t>& match_pos, uint_t match_length);

    // Converts the rare match map to pairs for further processing.
    RareMatchPairs convertMapToPairs(const RareMatchMap& rare_match_map);

    // Finds optimal pairs from given rare match pairs based on specific criteria.
    RareMatchPairs findOptimalPairs(const RareMatchPairs& rare_match_pairs);

    uint_t getMinMatchLength(std::vector<uint_t>& match_pos);

public:
    // Constructor initializes the finder with concatenated data and associated arrays.
    explicit RareMatchFinder(unsigned char* _concat_data, std::vector<uint_t>& _SA, std::vector<int_t>& _LCP, std::vector<int_da>& _DA, uint_t _first_seq_start, uint_t _first_seq_len, uint_t _second_seq_start, uint_t _second_seq_len);

    // Finds rare matches up to a specified maximum count.
    RareMatchPairs findRareMatch(uint_t max_match_count = 100);
};
