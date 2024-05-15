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
#include "rare_match.h"

// Function to save RareMatchPairs to a CSV file.
void saveRareMatchPairsToCSV(const RareMatchPairs& pairs, const std::string& filename, uint_t fst_len) {
    std::ofstream file(filename); // Opens the file for writing.
    if (!file.is_open()) { // Checks if the file is successfully opened.
        logger.error() << "Failed to open file: " << filename << std::endl; // Logs error if file cannot be opened.
        return; // Exits the function if file cannot be opened.
    }

    // Writes the header row of the CSV file.
    file << "Index,FirstPos,SecondPos,MatchLength,Weight\n";

    // Iterates through each RareMatchPair and writes its data to the CSV file.
    for (size_t i = 0; i < pairs.size(); ++i) {
        const auto& pair = pairs[i]; // References the current pair.
        // Writes the index and details of the current pair to the file.
        file << i + 1 << ","
            << pair.first_pos << ","
            << pair.second_pos - fst_len - 1 << ","
            << pair.match_length << ","
            << pair.weight << "\n";
    }

    file.close(); // Closes the file after writing.
    logger.info() << filename << " has been saved." << std::endl; // Logs successful save.
}

// Function to read RareMatchPairs from a CSV file, adjusting the second position based on a given offset (fst_len).
RareMatchPairs readRareMatchPairsFromCSV(const std::string& filename, uint_t fst_len) {
    std::ifstream file(filename); // Attempt to open the specified file for reading.
    RareMatchPairs pairs; // Vector to hold the read RareMatchPairs.

    // Check if the file was successfully opened.
    if (!file.is_open()) {
        logger.error() << "Failed to open file: " << filename << std::endl; // Log error message if file opening fails.
        return pairs; // Return the empty vector if the file cannot be opened.
    }

    std::string line;
    std::getline(file, line); // Read and discard the header line.

    // Read lines from the file one by one.
    while (std::getline(file, line)) {
        std::stringstream ss(line); // Use stringstream to parse the line.
        RareMatchPair pair; // Temporary object to hold the parsed data.
        std::string temp; // Temporary string to hold the value read from the stringstream.
        size_t i = 0; // Counter to track the current field being read.

        // Parse the current line by commas.
        while (std::getline(ss, temp, ',')) {
            switch (i++) {
            case 0: break; // Skip the index column.
            case 1: pair.first_pos = std::stoull(temp); break; // Parse the first position.
            case 2: pair.second_pos = std::stoull(temp) + fst_len + 1; break; // Parse and adjust the second position.
            case 3: pair.match_length = std::stoull(temp); break; // Parse the match length.
            case 4: pair.weight = std::stoull(temp); break; // Parse the weight.
            }
        }

        pairs.push_back(pair); // Add the parsed pair to the vector.
    }

    file.close(); // Close the file.
    return pairs; // Return the vector of parsed pairs.
}




LCPInterval::LCPInterval(const std::vector<int_t>& LCP_array, uint_t interval_size)
    : LCP(LCP_array), interval_size(interval_size), left(0), right(interval_size - 1), min_LCP_value(U_MAX) {
    if (interval_size == 1 && !LCP.empty()) {
        // Directly assign the value when interval size is 1
        min_LCP_value = LCP[0];
    }
    else {
        // Use deque to find the minimum LCP value for intervals larger than 1
        for (uint_t i = 0; i < interval_size && i < LCP.size(); ++i) {
            while (!min_deque.empty() && LCP[i] < LCP[min_deque.back()]) {
                min_deque.pop_back();
            }
            min_deque.push_back(i);
        }
        if (!min_deque.empty()) {
            min_LCP_value = LCP[min_deque.front()];
        }
    }
}

void LCPInterval::slideRight() {
    if (right + 1 < LCP.size()) {
        ++left;
        ++right;

        if (interval_size == 1) {
            // For interval size of 1, directly update min_LCP_value without using deque
            min_LCP_value = LCP[right];
        }
        else {
            // Update the deque for the new interval when interval size is greater than 1
            while (!min_deque.empty() && min_deque.front() < left) {
                min_deque.pop_front();
            }
            while (!min_deque.empty() && LCP[right] < LCP[min_deque.back()]) {
                min_deque.pop_back();
            }
            min_deque.push_back(right);

            // Update the min_LCP_value
            if (!min_deque.empty()) {
                min_LCP_value = LCP[min_deque.front()];
            }
        }
    }
}


// Returns the minimum LCP value within the current interval.
int_t LCPInterval::getMinLCP() const {
    return min_LCP_value;
}

// Determines if the current interval qualifies as a 'rare interval' based on specific criteria.
bool LCPInterval::isRareInterval() const {
    // Interval is not rare if the LCP value right before the start is greater or equal to the minimum LCP value...
    if (left > 0 && static_cast<uint_t>(LCP[left - 1]) >= min_LCP_value)
        return false;
    // ...or if the LCP value right after the end is greater or equal to the minimum LCP value.
    if (right < static_cast<uint_t>(LCP.size()) - 1 && static_cast<uint_t>(LCP[right + 1]) >= min_LCP_value)
        return false;
    return true;
}

// Checks if the right boundary of the interval has reached the end of the LCP array.
bool LCPInterval::isRightAtEnd() const {
    return right == LCP.size() - 1;
}

// Retrieves the left and right boundaries of the current interval.
std::pair<uint_t, uint_t> LCPInterval::getboundary() const {
    return std::make_pair(left, right);
}


// Constructor for the RareMatchFinder class.
// Initializes the class members with provided parameters and calculates additional properties.
RareMatchFinder::RareMatchFinder(unsigned char* _concat_data, // Pointer to the concatenated data array
    std::vector<uint_t>& _SA, // Suffix Array
    std::vector<int_t>& _LCP, // Longest Common Prefix array
    std::vector<int_da>& _DA, // Document Array indicating which sequence a suffix belongs to
    uint_t _first_seq_start, // Start position of the first sequence in the concatenated data
    uint_t _first_seq_len, // Length of the first sequence
    uint_t _second_seq_start, // Start position of the second sequence in the concatenated data
    uint_t _second_seq_len) // Length of the second sequence
    : concat_data(_concat_data),
    SA(_SA),
    LCP(_LCP),
    DA(_DA),
    first_seq_start(_first_seq_start),
    first_seq_len(_first_seq_len),
    second_seq_start(_second_seq_start),
    second_seq_len(_second_seq_len) {
    // Compute the minimum sequence length between the first and second sequences.
    min_seq_len = getMinValue(first_seq_len, second_seq_len);
    // Calculate the total length of the concatenated sequence.
    concat_seq_len = SA.size();
}

// Function to determine the minimum match length from the positions of matches.
// Takes a vector of match positions and returns the smallest match length.
uint_t RareMatchFinder::getMinMatchLength(std::vector<uint_t>& match_pos) {
    uint_t min_match_length = U_MAX; // Initialize with maximum possible value.
    for (auto& pos : match_pos) { // Iterate through each match position.
        if (pos >= second_seq_start) {
            // If the match position is within the second sequence, calculate the match length accordingly.
            min_match_length = getMinValue(min_match_length, second_seq_start + second_seq_len - pos);
        }
        else {
            // If the match position is within the first sequence, calculate the match length accordingly.
            min_match_length = getMinValue(min_match_length, first_seq_start + first_seq_len - pos);
        }
    }
    return min_match_length; // Return the smallest match length found.
}


// Finds rare matches up to a specified maximum count within the LCP array.
RareMatchPairs RareMatchFinder::findRareMatch(uint_t max_match_count) {
    // Limit the maximum match count to the minimum sequence length for efficiency.
    max_match_count = getMinValue(max_match_count, min_seq_len);
    uint_t lcp_interval_size = 0;
    bool is_match_found = false;
    RareMatchMap rare_match_map; // Stores unique rare matches.

    // Iterate through increasing LCP interval sizes until a match is found or max count is reached.
    while (!is_match_found && lcp_interval_size < max_match_count && ++lcp_interval_size) {
        LCPInterval lcp_interval(LCP, lcp_interval_size); // Create an LCP interval.

        // Process the LCP interval to find rare matches.
        while (!lcp_interval.isRightAtEnd()) {
            if (lcp_interval.isRareInterval()) { // Check for rare interval condition.
         
                std::pair<uint_t, uint_t> boundary = lcp_interval.getboundary();

                // Get match positions and types.
                std::vector<uint_t> match_pos;
                std::vector<bool> pos_type;
                getMatchPosAndType(boundary, match_pos, pos_type);

                // Calculate match length and get boundary for the current interval.
                uint_t match_length = getMinValue((uint_t)lcp_interval.getMinLCP(), getMinMatchLength(match_pos));

                // Convert match positions to RareMatch object and check counts.
                RareMatch rare_match(match_length, match_pos, pos_type);
                if (rare_match.first_count > 0 && rare_match.second_count > 0) {
                    is_match_found = true; // A match is found.
                    uint_t min_position = rare_match.min_key;
                    // Update or add the rare match in the map.
                    auto it = rare_match_map.find(min_position);
                    if (it != rare_match_map.end()) {
                        if(it->second.match_length < match_length)
                            it->second = rare_match;
                    } else {
                        rare_match_map[min_position] = rare_match;
                    }
                }
            }
            lcp_interval.slideRight(); // Move to the next interval.
        }
    }

    // Expand rare matches to the left and convert them to pairs.
    leftExpandRareMatchMap(rare_match_map);
    RareMatchPairs rare_match_pairs = convertMapToPairs(rare_match_map);

    // Sort and find optimal pairs from rare matches.
    std::sort(rare_match_pairs.begin(), rare_match_pairs.end());
    RareMatchPairs optimal_pairs = findOptimalPairs(rare_match_pairs);
    return optimal_pairs; // Return the optimal pairs found.
}


// Fills match positions and their types based on the specified boundary within the suffix array.
void RareMatchFinder::getMatchPosAndType(std::pair<uint_t, uint_t> boundary, // The boundary within the suffix array to consider
    std::vector<uint_t>& match_pos, // Vector to be filled with match positions
    std::vector<bool>& pos_type) { // Vector to be filled with position types (indicating sequence origin)
    uint_t left = boundary.first; // Starting position in the suffix array
    uint_t right = boundary.second; // Ending position in the suffix array

    // Adjust the left boundary to include the preceding element, if any
    if (left > 0)
        left--;

    // Iterate through the specified range in the suffix array
    for (uint_t i = left; i <= right; i++) {
        // Add the suffix array position to match positions
        match_pos.emplace_back(SA[i]);
        // Add the document array position to indicate the match's origin sequence
        pos_type.emplace_back(DA[i]);
    }

    return; // Explicit return for clarity
}

// Expands the matches in the rare match map towards the left to increase the match lengths.
void RareMatchFinder::leftExpandRareMatchMap(RareMatchMap& rare_match_map) { // A map of rare matches to be expanded
    // Iterate over each rare match pair in the map
    for (auto& pair : rare_match_map) {
        RareMatch& rare_match = pair.second; // Reference to the current rare match
        // Increase the match's length by expanding it towards the left
        // leftExpand returns the number of positions by which the match can be expanded
        rare_match.match_length += leftExpand(rare_match.match_pos, rare_match.match_length);
    }
}


uint_t RareMatchFinder::leftExpand(std::vector<uint_t>& match_pos, uint_t match_length) {
    if (match_pos.empty()) return 0; // Guard against empty input.

    uint_t max_expand_length = U_MAX;

    for (auto& pos : match_pos) {
        if (pos >= second_seq_start) {
            max_expand_length = getMinValue(max_expand_length, pos - second_seq_start);
        }
        else {
            max_expand_length = getMinValue(max_expand_length, pos - first_seq_start);
        } 
    }

    uint_t expand_length = 0;
    bool all_char_same = true;

    while (expand_length < max_expand_length) {
        expand_length++; // Move increment here to ensure correct initial value usage.
        if (match_pos[0] < expand_length) {
            all_char_same = false;
            expand_length--;
            break; // Ensure we don't go before start of concat_data.
        }
        unsigned char cur_char = concat_data[match_pos[0] - expand_length];
        for (uint_t i = 1; i < match_pos.size(); ++i) {
            // Check bounds for each position.
            if (match_pos[i] < expand_length || concat_data[match_pos[i] - expand_length] != cur_char) {
                all_char_same = false;
                break;
            }
        }

        if (!all_char_same) {
            expand_length--; // Undo the last increment if characters are not the same.
            break;
        }
    }

    // Update match_pos elements correctly using reference.
    for (uint_t& pos : match_pos) {
        pos -= expand_length; // Ensure this operation does not cause underflow.
    }

    return expand_length;
}


// Converts the map of rare matches into pairs for further processing.
// Each pair consists of positions from the first and second sequences, the match length, and a calculated weight.
RareMatchPairs RareMatchFinder::convertMapToPairs(const RareMatchMap& rare_match_map) {
    RareMatchPairs pairs;

    // Iterate through each rare match entry in the map.
    for (const auto& entry : rare_match_map) {
        const RareMatch& match = entry.second;

        // Separate positions based on whether they belong to the first or second sequence.
        std::vector<uint_t> first_seq_positions;
        std::vector<uint_t> second_seq_positions;
        for (uint_t i = 0; i < match.match_pos.size(); ++i) {
            if (!match.pos_type[i]) { // If it's from the first sequence
                first_seq_positions.emplace_back(match.match_pos[i]);
            }
            else { // From the second sequence
                second_seq_positions.emplace_back(match.match_pos[i]);
            }
        }

        // Create pairs from every combination of positions from the first and second sequences.
        for (auto first_pos : first_seq_positions) {
            for (auto second_pos : second_seq_positions) {
                // The weight of each pair is calculated based on the match length and the counts of matches in each sequence.
                // double weight = match.match_length / (match.first_count * match.second_count);
                double weight = match.match_length / getMinValue(match.first_count, match.second_count);
                /*if ((first_pos - this->first_seq_start) > (second_pos - this->second_seq_start)) {
                    weight -= ((first_pos - this->first_seq_start) - (second_pos - this->second_seq_start)) / (this->first_seq_len + this->second_seq_len);
                } else {
                    weight -= ((second_pos - this->second_seq_start) - (first_pos - this->first_seq_start)) / (this->first_seq_len + this->second_seq_len);
                }
                weight = getMaxValue(weight, 0.0);*/
                pairs.emplace_back(RareMatchPair{
                    first_pos,
                    second_pos,
                    match.match_length,
                    weight
                    });
            }
        }
    }

    return pairs; // Return the resulting vector of rare match pairs.
}


double gapCost(uint_t first_gap, uint_t second_gap) {
    if (first_gap - second_gap == 0)
        return 0;
    else
        return 2 * std::log2(std::abs(static_cast<int>(first_gap - second_gap)) + 1);
        /*return 0.01 * avg_seed_length * std::abs(static_cast<int>(first_gap - second_gap)) + 0.5 * std::log2(std::abs(static_cast<int>(first_gap - second_gap)) + 1);*/
}

// Finds the optimal sequence of rare match pairs based on a dynamic programming approach that maximizes the total weight.
RareMatchPairs RareMatchFinder::findOptimalPairs(const RareMatchPairs& rare_match_pairs) {
    if (rare_match_pairs.empty()) {
        return {}; // Return an empty vector if input is empty.
    }

    // Initialize vectors for scores and backtracks for dynamic programming.
    std::vector<double> scores(rare_match_pairs.size(), 0);
    std::vector<int_t> backtracks(rare_match_pairs.size(), -1);
    scores[0] = rare_match_pairs[0].weight; // Base case: first pair's score is its weight.

    // Iterate over each pair to calculate scores and find backtracks.
    for (uint_t i = 1; i < rare_match_pairs.size(); ++i) {
        scores[i] = rare_match_pairs[i].weight; // Start with the current pair's weight.
        for (int_t j = i - 1; j >= 0; --j) {
            
            // Check if current pair can follow the pair at j without overlap and if it improves the score.
            if (rare_match_pairs[i].first_pos >= rare_match_pairs[j].first_pos + rare_match_pairs[j].match_length &&
                rare_match_pairs[i].second_pos >= rare_match_pairs[j].second_pos + rare_match_pairs[j].match_length) {
                uint_t first_gap = rare_match_pairs[i].first_pos - rare_match_pairs[j].first_pos;
                uint_t second_gap = rare_match_pairs[i].second_pos - rare_match_pairs[j].second_pos;
                double gap_cost = gapCost(first_gap, second_gap);
                double new_score = scores[j];
                if(rare_match_pairs[i].weight - gap_cost > 0)
                    new_score = scores[j] + rare_match_pairs[i].weight - gap_cost;
                else
                    new_score = scores[j] + 0.1;
                if (new_score > scores[i]) {
                    scores[i] = new_score; // Update score if it's higher with the current combination.
                    backtracks[i] = j; // Record the backtrack index.
                }
            }
        }
    }

    // Find the index of the maximum score to start backtracking from.
    uint_t max_index = std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));

    // Reconstruct the optimal sequence of pairs using backtracks.
    RareMatchPairs optimal_pairs;
    for (int_t i = max_index; i != -1; i = backtracks[i]) {
        optimal_pairs.emplace_back(rare_match_pairs[i]);
    }

    // Reverse the sequence to get the correct order.
    std::reverse(optimal_pairs.begin(), optimal_pairs.end());
    return optimal_pairs; // Return the optimal sequence of pairs.
}
