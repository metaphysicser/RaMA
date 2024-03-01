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

#include "common.h"
#include "logging.h"
#include "utils.h"
#include "gsacak.h"
#include "rare_match.h"
#include "threadpool.h"
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include <algorithm>
#if defined(__GNUC__) || defined(__clang__)
#define CTZ(x) __builtin_ctz(x)
#elif defined(_MSC_VER)
static __inline unsigned long CTZ(unsigned long mask) {
    unsigned long index;
    _BitScanForward(&index, mask);
    return index;
}
#else
#error "Compiler not supports CTZ function!"
#endif
#define MAXM 20

#define ANCHORFINDER_NAME "anchorfinder.bin"

using Interval = std::pair<uint_t, uint_t>; // <start_pos, length>
using Intervals = std::vector<std::pair<Interval, Interval>>;

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

        for (Anchor* child : children) {
            // Recursively merge child pairs
            RareMatchPairs childPairs = child->mergeRareMatchPairs();
            mergedPairs.insert(mergedPairs.end(), childPairs.begin(), childPairs.end());
        }

        // Directly append rare_match_pairs if needed
        mergedPairs.insert(mergedPairs.end(), rare_match_pairs.begin(), rare_match_pairs.end());

        return mergedPairs;
    }

    // Destructor deletes all child Anchors
    ~Anchor() {
        for (Anchor* child : children) {
            delete child;
        }
    }
};


// Range Min Query
class RMQ : public Serializable {
private:
    uint_t N, block_size, block_num;
    int_t* LCP;
    std::vector<std::vector<uint_t>> st; // Sparse table
    std::vector<uint_t> pow, log; // Power and logarithm tables for fast computations
    std::vector<uint_t> pre, sub; // Precomputed values for block and sub-block queries
    std::vector<uint_t> belong, pos; // Auxiliary vectors for block decomposition
    std::vector<uint_t> f; // Auxiliary vector for queries

    // Builds the sparse table for RMQ
    void buildST();

    // Builds the precomputed tables for block and sub-block queries
    void buildSubPre();
    void buildSubPreParallel(); // Parallel version

    // Builds blocks for the RMQ structure
    void buildBlock();
    void buildBlockParallel(); // Parallel version

    // Utility functions for block decomposition
    int_t getBelong(int_t i) const;
    int_t getPos(int_t i) const;

public:
    // Default constructor initializes members
    explicit RMQ() : N(0), block_size(0), block_num(0), LCP(nullptr) {}

    // Constructor initializes LCP array and builds RMQ structure
    explicit RMQ(int_t* A, uint_t n, bool use_parallel = true);

    void setLCP(int_t* A);

    // Queries the minimum value in the range [l, r]
    int_t queryMin(uint_t l, uint_t r) const;

    // Serializes the RMQ structure to an output stream
    void serialize(std::ostream& out) const override;

    // Deserializes the RMQ structure from an input stream
    void deserialize(std::istream& in) override;
};


class AnchorFinder : public Serializable {
private:
    bool use_parallel; // Indicates whether to use parallel processing

    unsigned char* concat_data; // Concatenated sequence data

    uint_t concat_data_length; // Total length of the concatenated data
    uint_t first_seq_len; // Length of the first sequence
    uint_t second_seq_len; // Length of the second sequence

    uint_t* SA; // Suffix Array

    int_t* LCP; // Longest Common Prefix array

    int_da* DA; // Document Array indicating the origin sequence of each suffix

    uint_t* ISA; // Inverse Suffix Array

    RMQ rmq; // Range Minimum Query structure for LCP queries

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
    void constructISAParallel();

    // Locates anchors using a given thread pool, recursive depth, and intervals
    void locateAnchor(ThreadPool& pool, uint_t depth, uint_t task_id, Anchor* root, Interval first_interval, Interval second_interval);

    // Converts RareMatchPairs to Intervals considering specified intervals
    static Intervals RareMatchPairs2Intervals(const RareMatchPairs& rare_match_pairs, Interval first_interval, Interval second_interval, uint_t fst_length);

    // Converts a global index to a local index relative to concatenated sequences
    static uint_t indexFromGlogalToLocal(uint_t index, uint_t fst_length);

public:
    // Constructor initializes AnchorFinder with sequence data and optional parallel processing
    explicit AnchorFinder(std::vector<SequenceInfo>& data, bool use_parallel = true, std::string save_file_path = "./save/anchorfinder.bin", bool load_from_disk = false, bool save_to_disk = true);

    // Destructor cleans up allocated resources
    ~AnchorFinder();

    // Launches the anchor searching process
    RareMatchPairs lanuchAnchorSearching();
};
