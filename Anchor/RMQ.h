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
#include "config.h"
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
#define MAXM 32

// Range Min Query
class RMQ : public Serializable {
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


class SparseTable {
public:
    // Default constructor initializes members
    explicit SparseTable() : N(0) {}

    SparseTable(const int_t* LCP, size_t N) : N(N) {
        build(LCP);
    }

    int_t queryMin(size_t L, size_t R) const {
        if (L > R) std::swap(L, R);
        int_t j = log2[R - L + 1];
        return std::min(st[L][j], st[R - (1 << j) + 1][j]);
    }

private:
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

