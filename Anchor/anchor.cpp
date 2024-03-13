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

#include "anchor.h"

void saveIntervalsToCSV(const Intervals& intervals, const std::string& filename) {
    std::ofstream file(filename); // Opens the file for writing.
    if (!file.is_open()) { // Checks if the file is successfully opened.
        // Replace logger.error() with your logging mechanism or standard error output
        std::cerr << "Failed to open file: " << filename << std::endl; // Logs error if file cannot be opened.
        return; // Exits the function if file cannot be opened.
    }

    // Writes the header row of the CSV file.
    file << "Index,FirstStart,FirstLength,SecondStart,SecondLength\n";

    // Iterates through each pair in the Intervals and writes its data to the CSV file.
    for (size_t i = 0; i < intervals.size(); ++i) {
        const auto& intervalPair = intervals[i]; // References the current interval pair.
        const auto& firstInterval = intervalPair.first; // References the first interval in the pair.
        const auto& secondInterval = intervalPair.second; // References the second interval in the pair.

        // Writes the index and details of the current interval pair to the file.
        file << i + 1 << ","
            << firstInterval.first << "," // First interval's start position
            << firstInterval.second << "," // First interval's length
            << secondInterval.first << "," // Second interval's start position
            << secondInterval.second << "\n"; // Second interval's length
    }

    file.close(); // Closes the file after writing.
}

void RMQ::buildST() {
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


void RMQ::buildSubPre() {
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

void RMQ::buildSubPreParallel() {
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

// Sequentially constructs block information for RMQ
// This method implements a block-based approach to precompute RMQ information for each block.
// It utilizes a monotone stack to maintain the indices where LCP values are strictly decreasing.
// The f array is used to store block-wise precomputed RMQ data using bit manipulation.
void RMQ::buildBlock() {
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
//This method parallelizes the block-based RMQ preprocessing by dividing the sequence into chunks
//and processing each chunk in parallel, reducing overall computation time on multicore systems.
void RMQ::buildBlockParallel() {
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

RMQ::RMQ(int_t* a, uint_t n, bool use_parallel) {
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


int_t RMQ::queryMin(uint_t l, uint_t r) const {
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
        return LCP[l + CTZ(f[r] >> pos[l]) - 1]; // Directly query within the block
    }
}


// Returns the block ID to which the i-th element belongs
int_t RMQ::getBelong(int_t i) const {
    return (i - 1) / block_size + 1;
}

// Returns the position of the i-th element within its block
int_t RMQ::getPos(int_t i) const {
    return (i - 1) % block_size;
}


// Serializes the RMQ object state to an output stream.
// This includes saving the fundamental configurations and the various precomputed vectors
// necessary for quick RMQ queries.
void RMQ::serialize(std::ostream& out) const {
    // Save basic RMQ configuration numbers
    saveNumber(out, N);
    saveNumber(out, block_size);
    saveNumber(out, block_num);

    // Save precomputed vectors for RMQ algorithm
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


// Deserializes the RMQ object state from an input stream.
// This method loads the RMQ configuration and the precomputed vectors
// to fully reconstruct the RMQ state for future queries.
void RMQ::deserialize(std::istream& in) {
    // Load basic RMQ configuration numbers
    loadNumber(in, N);
    loadNumber(in, block_size);
    loadNumber(in, block_num);

    // Load precomputed vectors for RMQ algorithm
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

void RMQ::setLCP(int_t* A) {
    this->LCP = A;
}

// Constructor for AnchorFinder class
AnchorFinder::AnchorFinder(std::vector<SequenceInfo>& data, bool use_parallel, std::string save_file_path, bool load_from_disk, bool save_to_disk):
    use_parallel(use_parallel){
    first_seq_len = data[0].seq_len;
    second_seq_len = data[1].seq_len;
    concatSequence(data); // Concatenate sequences from input data
    logger.info() << "The concated data length is " << concat_data_length << std::endl;

    // Allocate memory for Suffix Array (SA), Longest Common Prefix (LCP), and Document Array (DA)
    this->SA = (uint_t*)malloc(concat_data_length * sizeof(uint_t));
    if (!SA) {
        logger.error() << "Failed to allocate " << concat_data_length * sizeof(uint_t) << "bytes of SA." << std::endl;
        exit(EXIT_FAILURE);
    }

    this->LCP = (int_t*)malloc(concat_data_length * sizeof(int_t));
    if (!LCP) {
        logger.error() << "Failed to allocate " << concat_data_length * sizeof(int_t) << "bytes of LCP." << std::endl;
        exit(EXIT_FAILURE);
    }

    this->DA = (int_da*)malloc(concat_data_length * sizeof(int_da));
    if (!DA) {
        logger.error() << "Failed to allocate " << concat_data_length * sizeof(int_t) << "bytes of DA." << std::endl;
        exit(EXIT_FAILURE);
    }

    this->ISA = (uint_t*)malloc(concat_data_length * sizeof(uint_t));
    if (!ISA) {
        logger.error() << "Failed to allocate " << concat_data_length * sizeof(uint_t) << "bytes of ISA." << std::endl;
        exit(EXIT_FAILURE);
    }


    ensureDirExists(save_file_path);
    std::string save_file_name = joinPaths(save_file_path, ANCHORFINDER_NAME);

    // Load arrays from disk if specified, otherwise construct the suffix array
    if (load_from_disk && fileExists(save_file_name) && loadFromFile(save_file_name)) {
        logger.info() << "AnchorFinder is loaded from " + save_file_name << std::endl;
    }
    else {
        if (load_from_disk) {
            logger.info() << "Fail to load " << save_file_name << ", start to construct arrays!" << std::endl;
        }
        logger.info() << "The suffix array is constructing..." << std::endl;
        gsacak(concat_data, SA, LCP, DA, concat_data_length);
        logger.info() << "The suffix array construction is finished!" << std::endl;

        logger.info() << "The sparse table is constructing..." << std::endl;
        this->rmq = RMQ(LCP, concat_data_length, use_parallel);
        logger.info() << "The sparse table construction is finished!" << std::endl;

        if (use_parallel)
            constructISAParallel();
        else
            constructISA(0, concat_data_length - 1);

        if (save_to_disk) {
            if(saveToFile(save_file_name))
                logger.info() << "AnchorFinder is saved into " + save_file_name << std::endl;
            else
                logger.info() << "Fail to save " + save_file_name << std::endl;
        }
            
    }

    if (logger.isDebugEnabled()) {
        printDebugInfo(SA, LCP, DA, concat_data_length);
    }
}

// Destructor for AnchorFinder class
AnchorFinder::~AnchorFinder() {
    // Free allocated memory
    if (concat_data) delete[] concat_data;
    if (SA) free(SA);
    if (LCP) free(LCP);
    if (DA) free(DA);
    if (ISA) free(ISA);
}

// Concatenates sequences and prepares data for suffix array construction
void AnchorFinder::concatSequence(std::vector<SequenceInfo>& data) {
    // Calculate total length of concatenated string
    uint_t total_length = 0;
    for (uint_t i = 0; i < data.size(); i++) {
        total_length += data[i].seq_len + 1; // Include separator
    }

    total_length++;  // Add 1 for the terminating 0

    // Allocate memory for concatenated string
    concat_data = new unsigned char[total_length];

    if (!concat_data) {
        logger.error() << "concat sequence could not allocate enough space" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Concatenate all strings with 1 as separator and terminate with 0
    uint_t index = 0;
    for (const auto& s : data) {
        std::copy(s.sequence.begin(), s.sequence.end(), concat_data + index);
        index += s.seq_len;
        concat_data[index] = 1; // Use 1 as separator
        index++;
    }
    concat_data[total_length - 1] = 0; // Set terminating 0

    concat_data_length = total_length;
}


// Serializes the state of the AnchorFinder object to an output stream.
// This includes all fundamental configurations, concatenated sequence data,
// various arrays such as Suffix Array (SA), Longest Common Prefix (LCP), 
// Document Array (DA), Inverse Suffix Array (ISA), and the associated RMQ structure.
void AnchorFinder::serialize(std::ostream& out) const {
    // Save basic configuration numbers
    saveNumber(out, concat_data_length);
    saveNumber(out, first_seq_len);
    saveNumber(out, second_seq_len);

    // Save concatenated sequence data and associated arrays
    saveArray(out, concat_data, concat_data_length);
    saveArray(out, SA, concat_data_length);
    saveArray(out, LCP, concat_data_length);
    saveArray(out, DA, concat_data_length);
    saveArray(out, ISA, concat_data_length);

    // Serialize the RMQ structure
    this->rmq.serialize(out);
}


// Deserializes the state of the AnchorFinder object from an input stream.
// This method loads the AnchorFinder configuration, concatenated sequence data,
// various arrays (SA, LCP, DA, ISA), and reconstructs the RMQ state for future queries.
void AnchorFinder::deserialize(std::istream& in) {
    // Load basic configuration numbers
    loadNumber(in, concat_data_length);
    loadNumber(in, first_seq_len);
    loadNumber(in, second_seq_len);

    // Load concatenated sequence data and associated arrays
    loadArray(in, concat_data, concat_data_length);
    loadArray(in, SA, concat_data_length);
    loadArray(in, LCP, concat_data_length);
    loadArray(in, DA, concat_data_length);
    loadArray(in, ISA, concat_data_length);

    // Deserialize the RMQ structure and set the LCP array for RMQ queries
    this->rmq.deserialize(in);
    this->rmq.setLCP(LCP);
}


// Prints the debug information for SA, LCP, DA, and ISA arrays.
// This includes indexing and the values within each array for debugging purposes.
void AnchorFinder::printDebugInfo(const uint_t* SA, const int_t* LCP, const int_da* DA, uint_t concat_data_length) {
    std::stringstream s;

    // Print indices for reference
    s << " index: ";
    for (uint_t i = 0; i < concat_data_length; ++i) {
        s << std::setw(6) << std::left << i << " ";
    }
    logger.debug() << s.str() << std::endl;
    s.str("");

    // Print Suffix Array (SA) values
    s << "    SA: ";
    for (uint_t i = 0; i < concat_data_length; ++i) {
        s << std::setw(6) << std::left << SA[i] << " ";
    }
    logger.debug() << s.str() << std::endl;
    s.str("");

    // Print Longest Common Prefix (LCP) values
    s << "   LCP: ";
    for (uint_t i = 0; i < concat_data_length; ++i) {
        s << std::setw(6) << std::left << LCP[i] << " ";
    }
    logger.debug() << s.str() << std::endl;
    s.str("");

    // Print Document Array (DA) values
    s << "    DA: ";
    for (uint_t i = 0; i < concat_data_length; ++i) {
        s << std::setw(6) << std::left << DA[i] << " ";
    }
    logger.debug() << s.str() << std::endl;
    s.str("");

    // Print Inverse Suffix Array (ISA) values
    s << "   ISA: ";
    for (uint_t i = 0; i < concat_data_length; ++i) {
        s << std::setw(6) << std::left << ISA[i] << " ";
    }
    logger.debug() << s.str() << std::endl;
    s.str("");
}

// Constructs the Inverse Suffix Array (ISA) for the given range.
// ISA array is built by marking the index of each element of SA in the ISA array,
// which is useful for various string processing algorithms.
void AnchorFinder::constructISA(uint_t start, uint_t end) {
    for (uint_t i = start; i <= getMinValue(end, concat_data_length - 1); i++) {
        ISA[SA[i]] = i;
    }
}


// Constructs the Inverse Suffix Array (ISA) in parallel.
// This method divides the task of building the ISA array into smaller chunks
// and processes each chunk in parallel using a thread pool, improving performance on multicore systems.
void AnchorFinder::constructISAParallel() {
    ThreadPool pool(std::thread::hardware_concurrency()); // Create a thread pool with a thread for each core
    const uint_t part_size = getMaxValue(std::ceil(static_cast<double>(concat_data_length) / std::thread::hardware_concurrency()), 1.0); // Calculate chunk size
    uint_t start = 0;
    uint_t end = part_size - 1;
    for (uint_t i = 0; i < concat_data_length; i += part_size) { // Divide work into chunks
        start = i;
        end = getMinValue(start + part_size - 1, concat_data_length - 1); // Ensure end does not exceed array bounds
        pool.enqueue([this, start, end]() { // Enqueue each chunk as a task
            this->constructISA(start, end);
            });
    }
    pool.waitAllTasksDone(); // Wait for all tasks to complete
}


// Initiates the process of searching for anchors in the sequences.
// Depending on the configuration, it either launches a parallel search using a thread pool
// or executes a single-threaded search.
RareMatchPairs AnchorFinder::lanuchAnchorSearching() {
    logger.info() << "Begin to search anchors" << std::endl;
    ThreadPool pool(std::thread::hardware_concurrency()); // Use thread pool for potential parallel execution
    uint_t depth = 0;
    Anchor* root = new Anchor(depth); // Create root anchor node
    Interval first_interval = std::make_pair(0, first_seq_len); // Define interval for the first sequence
    Interval second_interval = std::make_pair(0, second_seq_len); // Define interval for the second sequence
    uint_t task_id = 0;
    if (use_parallel) {
        pool.enqueue([this, &pool, depth, task_id, root, first_interval, second_interval]() { 
            this->locateAnchor(pool, depth, task_id, root, first_interval, second_interval);
            });
        pool.waitAllTasksDone();
    }
    else {
        locateAnchor(pool, depth, task_id, root, first_interval, second_interval); // Fallback to sequential search
    }
    RareMatchPairs first_anchors = root->rare_match_pairs;
    saveRareMatchPairsToCSV(first_anchors, "/mnt/f/code/vs_code/RaMA/output/first_anchor.csv", first_seq_len);

    RareMatchPairs final_anchors = root->mergeRareMatchPairs(); // Merge rare match pairs from the root anchor
    // RareMatchPairs final_anchors = verifyAnchors(root->mergeRareMatchPairs()); // Merge rare match pairs from the root anchor
    saveRareMatchPairsToCSV(final_anchors, "/mnt/f/code/vs_code/RaMA/output/final_anchor.csv", first_seq_len);

    delete root; // Clean up the root anchor
    logger.info() << "Finish searching anchors" << std::endl;

    return final_anchors;
    // return final_anchors;
}

// Launches the process of locating anchors within given intervals of two sequences.
// The method explores the given intervals, constructs new arrays based on the ISA,
// sorts them, and finds rare matches to determine new intervals for further exploration.
void AnchorFinder::locateAnchor(ThreadPool& pool, uint_t depth, uint_t task_id, Anchor* root, Interval first_interval, Interval second_interval) {
    // Log the start of a new task with its depth and task ID for debugging.
    logger.debug() << "Task " << task_id << " of depth " << depth << " begins" << std::endl;

    // Calculate new depth for recursive calls.
    uint_t new_depth = depth + 1;

    // Extract starting points and lengths from the intervals.
    uint_t first_seq_start = first_interval.first;
    uint_t fst_len = first_interval.second;
    uint_t second_seq_start = second_interval.first + first_seq_len + 1;
    uint_t scd_len = second_interval.second;

    // Return early if either sequence segment is empty.
    if (fst_len == 0 || scd_len == 0) {
        return;
    }
        
    uint_t new_array_len = fst_len + scd_len;

    // Prepare arrays to hold new SA, LCP, and DA values.
    std::vector<uint_t> new_index_of_SA;
    new_index_of_SA.reserve(new_array_len);

    for (uint_t i = first_seq_start; i < first_seq_start + fst_len; i++) {
        new_index_of_SA.emplace_back(ISA[i]);
    }
    for (uint_t i = second_seq_start; i < second_seq_start + scd_len; i++) {
        new_index_of_SA.emplace_back(ISA[i]);
    }

    // Sort the new SA indices to maintain the order.
    std::sort(new_index_of_SA.begin(), new_index_of_SA.end());

    // Create and populate new SA, LCP, and DA arrays based on the sorted indices.
    std::vector<uint_t> new_SA;
    std::vector<int_t> new_LCP;
    std::vector<int_da> new_DA;
    new_SA.reserve(new_array_len);
    new_LCP.reserve(new_array_len);
    new_DA.reserve(new_array_len);

    if (!new_index_of_SA.empty()) {
        uint_t last_index = new_index_of_SA[0];
        new_SA.emplace_back(SA[last_index]);
        new_DA.emplace_back(DA[last_index]);
        new_LCP.emplace_back(0);
        
        for (size_t i = 1; i < new_index_of_SA.size(); ++i) {
            auto index = new_index_of_SA[i];
            new_SA.emplace_back(SA[index]);
            new_DA.emplace_back(DA[index]);
            new_LCP.emplace_back(rmq.queryMin(last_index + 1, index));
            last_index = index;
        }
    }

    // Initialize RareMatchFinder and find optimal rare match pairs.
    RareMatchFinder rare_match_finder(concat_data, new_SA, new_LCP, new_DA, first_seq_start, fst_len, second_seq_start,scd_len);
    RareMatchPairs optimal_pairs = rare_match_finder.findRareMatch(100);

    // Convert rare match pairs to intervals for further exploration.
    Intervals rare_match_intervals = RareMatchPairs2Intervals(optimal_pairs, first_interval, second_interval, this->first_seq_len);

    // Update the anchor's rare match pairs with the optimal ones found.
    root->rare_match_pairs = optimal_pairs;

    // Recursively explore further intervals with new anchors.
    uint_t new_task_id = 0;
    for (const auto& interval_pair : rare_match_intervals) {
        Interval new_first_interval = interval_pair.first;
        Interval new_second_interval = interval_pair.second;
        Anchor* new_anchor = new Anchor(new_depth, root);
        root->children.emplace_back(new_anchor);
        // Parallel or sequential execution based on configuration.
        if (use_parallel) {
            pool.enqueue([this, &pool, new_depth, new_task_id, new_anchor, new_first_interval, new_second_interval]() {
            this->locateAnchor(pool, new_depth, new_task_id, new_anchor, new_first_interval, new_second_interval);
                });
        }
        else {
            locateAnchor(pool, new_depth, new_task_id, new_anchor, new_first_interval, new_second_interval);
        }   
        new_task_id++;
    }
    // Log the end of the current task.
    logger.debug() << "Task " << task_id << " of depth " << depth << " ends" << std::endl;
    return;
}

// Converts rare match pairs to intervals for anchor finding. The function determines
// intervals between rare matches for further analysis.
Intervals AnchorFinder::RareMatchPairs2Intervals(const RareMatchPairs& rare_match_pairs, Interval first_interval, Interval second_interval, uint_t fst_length) {
    // Return early if there are no rare match pairs.
    if (rare_match_pairs.empty())
        return {};
    Intervals intervals; // Store the resulting intervals.

    // Initialize starting points and lengths for the first and second sequence segments.
    uint_t start1 = first_interval.first;
    uint_t seq1_length = first_interval.second;
    uint_t start2 = second_interval.first + fst_length + 1;
    uint_t seq2_length = second_interval.second;
    uint_t i = 0;
    // Iterate over each rare match pair to calculate intervals.
    for (const auto& pair : rare_match_pairs) {
        // Calculate the start and end positions for each match in both sequences.
        uint_t match_start1 = pair.first_pos;
        uint_t match_start2 = pair.second_pos;
        uint_t match_end1 = match_start1 + pair.match_length - 1;
        uint_t match_end2 = match_start2 + pair.match_length - 1;

        // If current start points are less than or equal to match start points,
        // add the interval to the list and update the start points.
        if (start1 <= match_start1 && start2 <= match_start2) {
            intervals.push_back({ {start1, match_start1 - start1}, {indexFromGlogalToLocal(start2, fst_length), match_start2 - start2} });
        }
        else {
            logger.error() << "There is conflict in final anchors" << std::endl;
        }

        start1 = match_end1 + 1;
        start2 = match_end2 + 1;
    }

    // Calculate and add the final interval after the last match.
    Interval end1, end2;
    if (start1 >= seq1_length) {
        start1--;
        end1 = { start1, 0 };
    }
    else {
        end1 = { start1, seq1_length - start1 };
    }
        
    if (start2 >= fst_length + 1 + seq2_length) {
        start2--;
        end2 = { indexFromGlogalToLocal(start2, fst_length), 0 };
    }
    else {
        end2 = { indexFromGlogalToLocal(start2, fst_length), fst_length + 1 + seq2_length - start2 };
    }
        
    intervals.push_back({ end1, end2 });
    
    return intervals;
}

// Converts global index to local index, adjusting for the concatenation offset.
uint_t AnchorFinder::indexFromGlogalToLocal(uint_t index, uint_t fst_length) {
    // Adjust index if it belongs to the second sequence.
    return index > fst_length ? index - fst_length - 1 : index;
}

RareMatchPairs AnchorFinder::verifyAnchors(const RareMatchPairs& rare_match_pairs) {
    if (rare_match_pairs.empty()) return rare_match_pairs;

    // Sort the input vector (if not already sorted)
    RareMatchPairs sorted_pairs(rare_match_pairs);
    std::sort(sorted_pairs.begin(), sorted_pairs.end());

    RareMatchPairs verified;
    RareMatchPair current = sorted_pairs[0];

    for (size_t i = 1; i < sorted_pairs.size(); ++i) {
        if (current.hasOverlap(sorted_pairs[i])) {
            logger.error() << "Error: Overlapping RareMatchPairs detected.\n";
            exit(EXIT_FAILURE);
        }
        else if (current.isAdjacent(sorted_pairs[i])) {
            current.mergeWith(sorted_pairs[i]);
        }
        else {
            verified.push_back(current);
            current = sorted_pairs[i];
        }
    }

    // add the last processed RareMatchPair
    verified.push_back(current);

    return verified;
}
