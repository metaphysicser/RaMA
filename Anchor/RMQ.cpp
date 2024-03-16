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

#include "RMQ.h"

void LinearSparseTable::buildST() {
    // Initialize the first level of the sparse table with minimum values within each block
    int_t cur = 0, id = 1;
    for (uint_t i = 1; i <= N; ++i) {
        st[id][0] = getMinValue(st[id][0], (uint_t)LCP[i - 1]);
        belong[i] = id;
        pos[i] = cur;
        if (++cur == static_cast<int_t>(block_size)) {
            cur = 0;
            ++id;
        }
    }
    // Build the rest of the sparse table for efficient range minimum queries
    for (uint_t i = 1; i <= log[block_num]; ++i) {
        for (int_t j = 1; j + pow[i] - 1 <= block_num; ++j) {
            st[j][i] = getMinValue(st[j][i - 1], st[j + pow[i - 1]][i - 1]);
        }
    }
}


void LinearSparseTable::buildSubPre() {
    // Precompute minimum values for LCP within each block
    for (uint_t i = 1; i <= N; ++i) {
        if (belong[i] != belong[i - 1])
            pre[i] = LCP[i - 1];
        else {
            pre[i] = getMinValue(pre[i - 1], (uint_t)LCP[i - 1]);
        }

    }
    // Precompute minimum values for LCP within each block in reverse order
    for (uint_t i = N; i >= 1; --i) {
        if (i + 1 > N || belong[i] != belong[i + 1])
            sub[i] = LCP[i - 1];
        else
            sub[i] = getMinValue(sub[i + 1], (uint_t)LCP[i - 1]);
    }
}

void LinearSparseTable::buildSubPreParallel() {
    // Parallel version of buildSubPre using a thread pool for concurrency
    ThreadPool pool(std::thread::hardware_concurrency());

    for (uint_t block = 0; block < block_num; ++block) {
        pool.enqueue([this, block]() {
            uint_t start = block * block_size + 1;
            uint_t end = std::min((block + 1) * block_size, N);

            for (uint_t i = start; i <= end; ++i) {
                if (i == start || belong[i] != belong[i - 1])
                    pre[i] = LCP[i - 1];
                else
                    pre[i] = getMinValue(pre[i - 1], (uint_t)LCP[i - 1]);
            }

            });
    }


    for (uint_t block = 0; block < block_num; ++block) {
        pool.enqueue([this, block]() {
            uint_t start = block * block_size + 1;
            uint_t end = std::min((block + 1) * block_size, N);

            for (int_t i = static_cast<int_t>(end); i >= static_cast<int_t>(start); --i) {
                if (i == static_cast<int_t>(end) || i + 1 > N || belong[i] != belong[i + 1])
                    sub[i] = LCP[i - 1];
                else
                    sub[i] = getMinValue(sub[i + 1], (uint_t)LCP[i - 1]);
            }
            });
    }

    pool.waitAllTasksDone();
}

// Sequentially constructs block information for LinearSparseTable
// This method implements a block-based approach to precompute LinearSparseTable information for each block.
// It utilizes a monotone stack to maintain the indices where LCP values are strictly decreasing.
// The f array is used to store block-wise precomputed LinearSparseTable data using bit manipulation.
void LinearSparseTable::buildBlock() {
    int_t top = 0; // Stack pointer
    std::vector<int_t> s(block_size + 1, 0); // Monotone stack
    uint64_t bit = 1;
    for (uint_t i = 1; i <= N; ++i) {
        // Reset stack for each new block
        if (pos[i] == 0) top = 0;
        else f[i] = f[i - 1];
        // Maintain monotonicity of the stack
        while (top > 0 && LCP[s[top] - 1] >= LCP[i - 1])
            f[i] &= ~(bit << pos[s[top--]]); // Use bit manipulation to update f
        s[++top] = i; // Push current index onto stack
        f[i] |= (bit << pos[i]); // Set bit corresponding to current position
    }
}

//Parallel version of buildBlock using a thread pool for concurrency
//This method parallelizes the block-based LinearSparseTable preprocessing by dividing the sequence into chunks
//and processing each chunk in parallel, reducing overall computation time on multicore systems.
void LinearSparseTable::buildBlockParallel() {
    ThreadPool pool(std::thread::hardware_concurrency()); // Create a thread pool

    for (uint_t start = 1; start <= N; start += block_size) {
        pool.enqueue([this, start] {
            uint_t end = getMinValue(start + block_size - 1, N); // Determine block end
            int_t top = 0; // Stack pointer
            std::vector<int_t> s(block_size + 1, 0); // Monotone stack
            uint64_t bit = 1;
            for (uint_t i = start; i <= end; ++i) {
                // Reset stack for each new block
                if (pos[i] == 0) top = 0;
                else f[i] = f[i - 1];
                // Maintain monotonicity of the stack
                while (top > 0 && LCP[s[top] - 1] >= LCP[i - 1])
                    f[i] &= ~(bit << pos[s[top--]]); // Use bit manipulation to update f
                s[++top] = i; // Push current index onto stack
                f[i] |= (bit << pos[i]); // Set bit corresponding to current position
            }
            });
    }

    pool.waitAllTasksDone(); // Wait for all tasks to complete
}

LinearSparseTable::LinearSparseTable(int_t* a, uint_t n, bool use_parallel) {
    LCP = a; // Directly use the provided array for LCP values
    N = n; // Set the total number of elements
    // Initialize vectors with appropriate sizes and default values
    belong.resize(N + 1, 0);
    pos.resize(N + 1, 0);
    pow.resize(MAXM, 0);
    log.resize(N + 1, 0);
    pre.resize(N + 1, 0);
    sub.resize(N + 1, 0);
    f.resize(N + 1, 0);

    // Calculate block size and number of blocks based on input size
    block_size = getMinValue((int)(log2(N) * 1.5), 63);
    block_num = (N + block_size - 1) / block_size;

    // Initialize 'pow' and 'log' arrays for fast range queries
    pow[0] = 1;
    for (uint_t i = 1; i < MAXM; ++i) pow[i] = pow[i - 1] * 2;
    for (uint_t i = 2; i <= block_num; ++i) log[i] = log[i / 2] + 1;

    // Initialize sparse table with maximum values
    st.resize(block_num + 1);
    for (uint_t i = 0; i < st.size(); ++i) {
        st[i].resize(log[block_num] + 1, U_MAX);
    }

    // Build the sparse table and preprocess LCP array
    buildST();
    if (use_parallel) { // Choose parallel or sequential preprocessing based on flag
        buildSubPreParallel();
        buildBlockParallel();
    }
    else {
        buildSubPre();
        buildBlock();
    }
}


int_t LinearSparseTable::queryMin(uint_t l, uint_t r) const {
    assert(l >= 0 && r <= N); // Ensure query indices are within bounds
    ++l; // Convert to 1-based indexing
    ++r;
    int_t bl = belong[l], br = belong[r]; // Get block IDs for l and r
    if (bl != br) { // If l and r are in different blocks
        int_t ans1 = I_MAX; // Initialize answer for block query
        if (br - bl > 1) { // Query for blocks between l and r
            int_t p = log[br - bl - 1];
            ans1 = getMinValue(st[bl + 1][p], st[br - pow[p]][p]);
        }
        int_t ans2 = getMinValue(sub[l], pre[r]); // Query for prefix and suffix within blocks
        return getMinValue(ans1, ans2); // Return the overall minimum
    }
    else { // If l and r are in the same block
        uint_t m = f[r] >> pos[l];
        uint_t q = CTZ(m);
        return LCP[l + CTZ(f[r] >> pos[l]) - 1]; // Directly query within the block
    }
}


// Returns the block ID to which the i-th element belongs
int_t LinearSparseTable::getBelong(int_t i) const {
    return (i - 1) / block_size + 1;
}

// Returns the position of the i-th element within its block
int_t LinearSparseTable::getPos(int_t i) const {
    return (i - 1) % block_size;
}


// Serializes the LinearSparseTable object state to an output stream.
// This includes saving the fundamental configurations and the various precomputed vectors
// necessary for quick LinearSparseTable queries.
void LinearSparseTable::serialize(std::ostream& out) const {
    // Save basic LinearSparseTable configuration numbers
    saveNumber(out, N);
    saveNumber(out, block_size);
    saveNumber(out, block_num);

    // Save precomputed vectors for LinearSparseTable algorithm
    saveVector(out, pow);
    saveVector(out, log);
    saveVector(out, pre);
    saveVector(out, sub);
    saveVector(out, belong);
    saveVector(out, pos);
    saveVector(out, f);

    // Save the 2D sparse table
    saveVector2D(out, st);
}


// Deserializes the LinearSparseTable object state from an input stream.
// This method loads the LinearSparseTable configuration and the precomputed vectors
// to fully reconstruct the LinearSparseTable state for future queries.
void LinearSparseTable::deserialize(std::istream& in) {
    // Load basic LinearSparseTable configuration numbers
    loadNumber(in, N);
    loadNumber(in, block_size);
    loadNumber(in, block_num);

    // Load precomputed vectors for LinearSparseTable algorithm
    loadVector(in, pow);
    loadVector(in, log);
    loadVector(in, pre);
    loadVector(in, sub);
    loadVector(in, belong);
    loadVector(in, pos);
    loadVector(in, f);

    // Load the 2D sparse table
    loadVector2D(in, st);
}

void LinearSparseTable::setLCP(int_t* A) {
    this->LCP = A;
}
