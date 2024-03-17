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
 // Created: 2024-03-14
#pragma once

#include "logging.h"
#include "utils.h"
#include "gsacak.h"
#include "rare_match.h"
#include "threadpool.h"
#include <algorithm>

#define MAXM 32

// An inline function to count trailing zeros (CTZ) in a 64-bit unsigned integer.
// The CTZ operation counts the number of zero bits starting from the least significant bit until the first one bit is found.
inline uint_t CTZ(uint64_t x) {
    // If the input number is 0, then technically all bits are trailing zeros.
    // Since we're working with a 64-bit number, return 64 to indicate this.
    if (x == 0) return 64;

    uint_t count = 0; // Initialize a counter to keep track of the number of trailing zeros.
    // Loop until we find a bit that is set to 1.
    while ((x & 1) == 0) {
        count++; // Increment the counter for each trailing zero found.
        x >>= 1; // Right-shift the number to check the next bit in the next iteration.
    }
    return count; // Return the total number of trailing zeros found.
}


class LinearSparseTable : public Serializable {
private:
    uint_t N, block_size, block_num;
    int_t* LCP;
    std::vector<std::vector<uint_t>> st; // Sparse table
    std::vector<uint_t> pow, log; // Power and logarithm tables for fast computations
    std::vector<uint_t> pre, sub; // Precomputed values for block and sub-block queries
    std::vector<uint_t> belong, pos; // Auxiliary vectors for block decomposition
    std::vector<uint64_t> f; // Auxiliary vector for queries

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
    explicit LinearSparseTable() : N(0), block_size(0), block_num(0), LCP(nullptr) {}

    // Constructor initializes LCP array and builds RMQ structure
    explicit LinearSparseTable(int_t* A, uint_t n, bool use_parallel = true);

    void setLCP(int_t* A);

    // Queries the minimum value in the range [l, r]
    int_t queryMin(uint_t l, uint_t r) const;

    // Serializes the RMQ structure to an output stream
    void serialize(std::ostream& out) const override;

    // Deserializes the RMQ structure from an input stream
    void deserialize(std::istream& in) override;
};


class SparseTable : public Serializable {
public:
    // Default constructor initializes members
    explicit SparseTable() : N(0) {}

    SparseTable(const int_t* LCP, size_t N) : LCP(LCP), N(N) {
        build(LCP);
    }

    int_t queryMin(size_t L, size_t R) const {
        if (L > R) std::swap(L, R);
        int_t j = log2[R - L + 1];
        return std::min(st[L][j], st[R - (1 << j) + 1][j]);
    }

    // Serializes the RMQ structure to an output stream
    void serialize(std::ostream& out) const override {
        saveNumber(out, N);
        saveVector2D(out, st);
        saveVector(out, log2);
    }

    // Deserializes the RMQ structure from an input stream
    void deserialize(std::istream& in) override {
        loadNumber(in, N);
        loadVector2D(in, st);
        loadVector(in, log2);
    }

    void setLCP(int_t* A) {
        this->LCP = A;
    }

private:
    const int_t* LCP;
    size_t N;
    std::vector<std::vector<int_t>> st;
    std::vector<int_t> log2;

    void build(const int_t* LCP) {
        size_t K = std::log2(N) + 1;
        st.assign(N, std::vector<int_t>(K, I_MAX));
        log2.assign(N + 1, 0);

        for (size_t i = 2; i <= N; i++) {
            log2[i] = log2[i / 2] + 1;
        }

        for (size_t i = 0; i < N; i++) {
            st[i][0] = LCP[i];
        }

        for (size_t j = 1; j <= K; j++) {
            for (size_t i = 0; i + (1 << j) <= N; i++) {
                st[i][j] = std::min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
            }
        }
    }


};

