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
 // Created: 2024-02-04
#pragma once

#include "gsacak.h"
#include "logging.h"
#include "utils.h"
#include "gsacak.h"
#include "rare_match.h"
#include "threadpool.h"
#include "RMQ.h"

#define SAVE_DIR "save"
#define ANCHORFINDER_NAME "anchorfinder.bin"
#define FIRST_ANCHOR_NAME "first_anchor.csv"
#define FINAL_ANCHOR_NAME "final_anchor.csv"

struct Interval
{
    uint_t pos1;
    uint_t len1;

    uint_t pos2;
    uint_t len2;

    Interval(uint_t p1, uint_t l1, uint_t p2, uint_t l2) : pos1(p1), len1(l1), pos2(p2), len2(l2) {}
    Interval() : pos1(0), len1(0), pos2(0), len2(0) {}
};

using Intervals = std::vector<Interval>;

void saveIntervalsToCSV(const Intervals& intervals, const std::string& filename); 

struct Anchor {
    uint_t depth;
    Anchor* parent;
    std::vector<Anchor*> children;
    RareMatchPairs rare_match_pairs;

    // Constructor initializes depth and parent
    Anchor(uint_t d = 0, Anchor* p = nullptr) : depth(d), parent(p) {}

    // Merges RareMatchPairs from this node and all descendants
    RareMatchPairs mergeRareMatchPairs() {
        RareMatchPairs mergedPairs;

        // Reserve memory to improve efficiency
        mergedPairs.reserve(100); 
        uint_t rare_pair_index = 0;

        for (Anchor* child : children) {
            // Recursively merge child pairs
            RareMatchPairs childPairs = child->mergeRareMatchPairs();
            mergedPairs.insert(mergedPairs.end(), childPairs.begin(), childPairs.end());
            if (rare_match_pairs.begin() + rare_pair_index != rare_match_pairs.end())
                mergedPairs.insert(mergedPairs.end(), rare_match_pairs.begin()+rare_pair_index, rare_match_pairs.begin() + rare_pair_index + 1);
            rare_pair_index++;
        }

        return mergedPairs;
    }

    // Destructor deletes all child Anchors
    ~Anchor() {
        for (Anchor* child : children) {
            delete child;
        }
    }
};



class AnchorFinder : public Serializable {
private:
    uint_t thread_num; // Indicates whether to use parallel processing

    std::string save_file_path; // Path to the directory where the AnchorFinder state is saved

    uint_t max_match_count; // Maximum number of rare matches to find

    unsigned char* concat_data; // Concatenated sequence data

    uint_t concat_data_length; // Total length of the concatenated data
    uint_t first_seq_len; // Length of the first sequence
    uint_t second_seq_len; // Length of the second sequence

    uint_t* SA; // Suffix Array

    int_t* LCP; // Longest Common Prefix array

    int_da* DA; // Document Array indicating the origin sequence of each suffix

    uint_t* ISA; // Inverse Suffix Array

    LinearSparseTable rmq; // Range Minimum Query structure for LCP queries

    // Concatenates sequences from provided data
    void concatSequence(std::vector<SequenceInfo>& data);

    // Prints debug information including SA, LCP, and DA
    void printDebugInfo(const uint_t* SA, const int_t* LCP, const int_da* DA, uint_t concat_data_length);

    // Serialization of AnchorFinder state to an output stream
    void serialize(std::ostream& out) const override;

    // Deserialization of AnchorFinder state from an input stream
    void deserialize(std::istream& in) override;

    // Constructs the Inverse Suffix Array (ISA) for a given range
    void constructISA(uint_t start, uint_t end);

    // Constructs the ISA in parallel
    void constructISAParallel(uint_t thread_num);

    // Locates anchors using a given thread pool, recursive depth, and intervals
    void locateAnchor(ThreadPool& pool, uint_t depth, uint_t task_id, Anchor* root, Interval interval);

    RareMatchPairs verifyAnchors(const RareMatchPairs& rare_match_pairs);


public:
    // Constructor initializes AnchorFinder with sequence data and optional parallel processing
    explicit AnchorFinder(std::vector<SequenceInfo>& data, std::string save_file_path, uint_t thread_num = 0, bool load_from_disk = false, bool save_to_disk = true, uint_t max_match_count = 100);

    // Destructor cleans up allocated resources
    ~AnchorFinder();

    // Launches the anchor searching process
    RareMatchPairs lanuchAnchorSearching();

    // Converts RareMatchPairs to Intervals considering specified intervals
    static Intervals RareMatchPairs2Intervals(const RareMatchPairs& rare_match_pairs, Interval interval, uint_t fst_length);

    // Converts a global index to a local index relative to concatenated sequences
    static uint_t indexFromGlogalToLocal(uint_t index, uint_t fst_length);
};
