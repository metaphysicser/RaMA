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

std::mutex mtx;  // 定义互斥锁
uint_t total_sub_suffix_array = 0;  // 初始化计数器

void increment_count(uint_t& total_sub_suffix_array, uint_t length) {
	std::lock_guard<std::mutex> lock(mtx);
	total_sub_suffix_array += length;
}


void saveIntervalsToCSV(const Intervals& intervals, const std::string& filename) {
	std::ofstream file(filename); // Opens the file for writing.
	if (!file.is_open()) { // Checks if the file is successfully opened.
		// Replace logger.error() with your logging mechanism or standard error output
		logger.error() << "Failed to open file: " << filename << std::endl; // Logs error if file cannot be opened.
		return; // Exits the function if file cannot be opened.
	}

	// Writes the header row of the CSV file.
	file << "Index,FirstStart,FirstLength,SecondStart,SecondLength\n";

	// Iterates through each pair in the Intervals and writes its data to the CSV file.
	for (size_t i = 0; i < intervals.size(); ++i) {
		const auto& interval = intervals[i]; // References the current interval pair.

		// Writes the index and details of the current interval pair to the file.
		file << i + 1 << ","
			<< interval.pos1 << "," // First interval's start position
			<< interval.len1 << "," // First interval's length
			<< interval.pos2 << "," // Second interval's start position
			<< interval.len2 << "\n"; // Second interval's length
	}

	file.close(); // Closes the file after writing.
	logger.info() << filename << " has been saved" << std::endl;
}

// Constructor for AnchorFinder class
AnchorFinder::AnchorFinder(std::vector<SequenceInfo>& data, std::string save_file_path, uint_t thread_num, bool load_from_disk, bool save_to_disk, uint_t max_match_count) :
	save_file_path(save_file_path),
	thread_num(thread_num),
	max_match_count(max_match_count) {
	first_seq_len = data[0].seq_len;
	second_seq_len = data[1].seq_len;
	concatSequence(data); // Concatenate sequences from input data
	logger.info() << "The concated data length is " << concat_data_length << std::endl;

	// Allocate memory for Suffix Array (SA), Longest Common Prefix (LCP), and Document Array (DA)
	this->SA = (uint_t*)malloc(concat_data_length * sizeof(uint_t));
	if (!SA) {
		logger.error() << "Failed to allocate " << concat_data_length * sizeof(uint_t) << "bytes of SA." << std::endl;
		logger.error() << "RaMA Exit!" << std::endl;
		exit(EXIT_FAILURE);
	}

	this->LCP = (int_t*)malloc(concat_data_length * sizeof(int_t));
	if (!LCP) {
		logger.error() << "Failed to allocate " << concat_data_length * sizeof(int_t) << "bytes of LCP." << std::endl;
		logger.error() << "RaMA Exit!" << std::endl;
		exit(EXIT_FAILURE);
	}

	this->DA = (int_da*)malloc(concat_data_length * sizeof(int_da));
	if (!DA) {
		logger.error() << "Failed to allocate " << concat_data_length * sizeof(int_t) << "bytes of DA." << std::endl;
		logger.error() << "RaMA Exit!" << std::endl;
		exit(EXIT_FAILURE);
	}

	this->ISA = (uint_t*)malloc(concat_data_length * sizeof(uint_t));
	if (!ISA) {
		logger.error() << "Failed to allocate " << concat_data_length * sizeof(uint_t) << "bytes of ISA." << std::endl;
		logger.error() << "RaMA Exit!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string bin_file_dir = joinPaths(save_file_path, SAVE_DIR);
	ensureDirExists(bin_file_dir);
	std::string save_file_name = joinPaths(bin_file_dir, ANCHORFINDER_NAME);

	// Load arrays from disk if specified, otherwise construct the suffix array
	if (load_from_disk && fileExists(save_file_name) && loadFromFile(save_file_name)) {
		logger.info() << "AnchorFinder is loaded from " + save_file_name << std::endl;
	}
	else {
		if (load_from_disk) {
			logger.info() << "Fail to load " << save_file_name << ", start to construct arrays!" << std::endl;
		}
		logger.info() << "The suffix array is constructing..." << std::endl;
		try {
			gsacak(concat_data, SA, LCP, DA, concat_data_length);
			logger.info() << "The suffix array construction is finished!" << std::endl;
		}
		catch (std::exception& e) {
			logger.error() << "Error: " << e.what() << std::endl;
			logger.error() << "RaMA Exit!" << std::endl;
			exit(EXIT_FAILURE);
		}

		logger.info() << "The sparse table is constructing..." << std::endl;
		try {
			this->rmq = LinearSparseTable(LCP, concat_data_length, thread_num);
			logger.info() << "The sparse table construction is finished!" << std::endl;
		}
		catch (std::exception& e) {
			logger.error() << "Error: " << e.what() << std::endl;
			logger.error() << "RaMA Exit!" << std::endl;
			exit(EXIT_FAILURE);
		}


		if (thread_num)
			constructISAParallel(thread_num);
		else
			constructISA(0, concat_data_length - 1);

		if (save_to_disk) {
			if (saveToFile(save_file_name))
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
void AnchorFinder::constructISAParallel(uint_t thread_num) {
	ThreadPool pool(thread_num); // Create a thread pool with a thread for each core
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
	total_sub_suffix_array = 0;
	ThreadPool pool(thread_num); // Use thread pool for potential parallel execution
	uint_t depth = 0;
	Anchor* root = new Anchor(depth); // Create root anchor node
	Interval interval(0, first_seq_len, 0, second_seq_len); // Define interval
	uint_t task_id = 0;
	if (thread_num) {
		pool.enqueue([this, &pool, depth, task_id, root, interval]() {
			this->locateAnchor(pool, depth, task_id, root, interval);
			});
		pool.waitAllTasksDone();
	}
	else {
		locateAnchor(pool, depth, task_id, root, interval); // Fallback to sequential search
	}
	RareMatchPairs first_anchors = root->rare_match_pairs;
	saveRareMatchPairsToCSV(first_anchors, joinPaths(save_file_path, FIRST_ANCHOR_NAME), first_seq_len);

	// RareMatchPairs final_anchors = root->mergeRareMatchPairs(); // Merge rare match pairs from the root anchor
	RareMatchPairs final_anchors = verifyAnchors(root->mergeRareMatchPairs()); // Merge rare match pairs from the root anchor
	saveRareMatchPairsToCSV(final_anchors, joinPaths(save_file_path, FINAL_ANCHOR_NAME), first_seq_len);

	logger.info() << "New sub suffix array length is " << total_sub_suffix_array - (first_seq_len + second_seq_len) << ". Compared to a multiple of the original sequence length is " << (float)(total_sub_suffix_array - (first_seq_len + second_seq_len)) / (first_seq_len + second_seq_len) << std::endl;
	delete root; // Clean up the root anchor
	logger.info() << "Finish searching anchors" << std::endl;

	return final_anchors;
}

// Launches the process of locating anchors within given intervals of two sequences.
// The method explores the given intervals, constructs new arrays based on the ISA,
// sorts them, and finds rare matches to determine new intervals for further exploration.
void AnchorFinder::locateAnchor(ThreadPool& pool, uint_t depth, uint_t task_id, Anchor* root, Interval interval) {
	// Log the start of a new task with its depth and task ID for debugging.
	logger.debug() << "Task " << task_id << " of depth " << depth << " begins" << std::endl;

	// Calculate new depth for recursive calls.
	uint_t new_depth = depth + 1;

	// Extract starting points and lengths from the intervals.
	uint_t first_seq_start = interval.pos1;
	uint_t fst_len = interval.len1;
	uint_t second_seq_start = interval.pos2 + first_seq_len + 1;
	uint_t scd_len = interval.len2;

	// Return early if either sequence segment is empty.
	if (fst_len == 0 || scd_len == 0) {
		return;
	}

	uint_t new_array_len = fst_len + scd_len;
	increment_count(total_sub_suffix_array, new_array_len);

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
	/*new_SA.reserve(new_array_len);
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
	}*/
	new_SA.resize(new_array_len);
	new_LCP.resize(new_array_len);
	new_DA.resize(new_array_len);

	if (!new_index_of_SA.empty()) {
		uint_t last_index = new_index_of_SA[0];
		new_SA[0] = SA[last_index];
		new_DA[0] = DA[last_index];
		new_LCP[0] = 0;
#pragma omp parallel for
		for (size_t i = 1; i < new_index_of_SA.size(); ++i) {
			auto index = new_index_of_SA[i];
			new_SA[i] = SA[index];
			new_DA[i] = DA[index];
			new_LCP[i] = rmq.queryMin(new_index_of_SA[i - 1] + 1, index);
		}

	}

	// Initialize RareMatchFinder and find optimal rare match pairs.
	RareMatchFinder rare_match_finder(concat_data, new_SA, new_LCP, new_DA, first_seq_start, fst_len, second_seq_start, scd_len);
	RareMatchPairs optimal_pairs = rare_match_finder.findRareMatch(max_match_count);

	if (optimal_pairs.empty())
		return;

	// Convert rare match pairs to intervals for further exploration.
	Intervals rare_match_intervals = RareMatchPairs2Intervals(optimal_pairs, interval, this->first_seq_len);

	// Update the anchor's rare match pairs with the optimal ones found.
	root->rare_match_pairs = optimal_pairs;

	// Recursively explore further intervals with new anchors.
	uint_t new_task_id = 0;

	for (const auto& new_interval : rare_match_intervals) {
		Anchor* new_anchor = new Anchor(new_depth, root);
		root->children.emplace_back(new_anchor);
		// Parallel or sequential execution based on configuration.
		if (thread_num) {
			pool.enqueue([this, &pool, new_depth, new_task_id, new_anchor, new_interval]() {
				this->locateAnchor(pool, new_depth, new_task_id, new_anchor, new_interval);
				});
		}
		else {
			locateAnchor(pool, new_depth, new_task_id, new_anchor, new_interval);
		}
		new_task_id++;
	}
	// Log the end of the current task.
	logger.debug() << "Task " << task_id << " of depth " << depth << " ends" << std::endl;
	return;
}

// Converts rare match pairs to intervals for anchor finding. The function determines
// intervals between rare matches for further analysis.
Intervals AnchorFinder::RareMatchPairs2Intervals(const RareMatchPairs& rare_match_pairs, Interval interval, uint_t fst_length) {
	Intervals intervals; // Store the resulting intervals.

	// Return early if there are no rare match pairs.
	if (rare_match_pairs.empty()) {
		intervals.emplace_back(interval);
		return intervals;
	}


	// Initialize starting points and lengths for the first and second sequence segments.
	uint_t start1 = interval.pos1;
	uint_t seq1_length = interval.len1;
	uint_t start2 = interval.pos2 + fst_length + 1;
	uint_t seq2_length = interval.len2;
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
			intervals.emplace_back(Interval(start1, match_start1 - start1, indexFromGlogalToLocal(start2, fst_length), match_start2 - start2));
		}
		else {
			logger.error() << "There is conflict in final anchors" << std::endl;
		}

		start1 = match_end1 + 1;
		start2 = match_end2 + 1;
	}

	// Calculate and add the final interval after the last match.
	Interval end;
	if (start1 >= seq1_length) {
		start1--;
		end.pos1 = start1;
		end.len1 = 0;
	}
	else {
		end.pos1 = start1;
		end.len1 = seq1_length - start1;
	}

	if (start2 >= fst_length + 1 + seq2_length) {
		start2--;
		end.pos2 = indexFromGlogalToLocal(start2, fst_length);
		end.len2 = 0;
	}
	else {
		end.pos2 = indexFromGlogalToLocal(start2, fst_length);
		end.len2 = fst_length + 1 + seq2_length - start2;
	}

	intervals.emplace_back(end);

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