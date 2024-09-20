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
 // Created: 2024-02-29

# include "pairwise_alignment.h"

// Converts a buffer of CIGAR operations (represented as compact integers) into a vector.
cigar convertToCigarVector(uint32_t* cigar_buffer, int cigar_length) {
	cigar result_cigar; // Initialize an empty vector to store the result.

	// Iterate through each element in the cigar_buffer.
	for (int i = 0; i < cigar_length; ++i) {
		// Add the current CIGAR operation from the buffer to the result vector.
		result_cigar.push_back(cigar_buffer[i]);
	}

	return result_cigar; // Return the populated vector of CIGAR operations.
}


uint32_t cigarToInt(char operation, uint32_t len) {
	uint32_t opCode;
	// Convert CIGAR operation character to an operation code
	switch (operation) {
	case 'M': opCode = 0x0; break; // Match
	case 'I': opCode = 0x1; break; // Insertion
	case 'D': opCode = 0x2; break; // Deletion
	case '=': opCode = 0x7; break; // Sequence match
	case 'X': opCode = 0x8; break; // Mismatch
		// Add cases for other SAM specification operation codes as needed
	default: opCode = 0xF; break; // Unknown operation
	}
	// Combine operation length and code into a single uint32_t value
	return (len << 4) | opCode; // Shift length left by 4 bits, then combine with opCode
}


void intToCigar(uint32_t cigar, char& operation, uint32_t& len) {
	uint32_t opCode = cigar & 0xF; // Extract the lower 4 bits as the operation code
	len = cigar >> 4; // Extract the length by shifting right by 4 bits

	// Convert operation code back to a CIGAR operation character
	switch (opCode) {
	case 0x0: operation = 'M'; break; // Match
	case 0x1: operation = 'I'; break; // Insertion
	case 0x2: operation = 'D'; break; // Deletion
	case 0x7: operation = '='; break; // Sequence match
	case 0x8: operation = 'X'; break; // Mismatch
		// Add cases for other SAM specification operation codes as needed
	default: operation = '?'; break; // Unknown operation
	}
}


PairAligner::PairAligner(std::string save_file_path, int_t match, int_t mismatch, int_t gap_open1, int_t gap_extension1, int_t gap_open2, int_t gap_extension2, uint_t thread_num) :
	save_file_path(save_file_path),
	match(match),
	mismatch(mismatch),
	gap_open1(gap_open1),
	gap_extension1(gap_extension1),
	gap_open2(gap_open2),
	gap_extension2(gap_extension2),
	thread_num(thread_num) {
	attributes = wavefront_aligner_attr_default;

	//attributes.distance_metric = gap_affine;
	//attributes.affine_penalties.mismatch = 2;      // X > 0
	//attributes.affine_penalties.gap_opening = 3;   // O >= 0
	//attributes.affine_penalties.gap_extension = 1; // E > 0

	attributes.distance_metric = gap_affine_2p;
	attributes.affine2p_penalties.match = match;
	attributes.affine2p_penalties.mismatch = mismatch;       // X > 0
	attributes.affine2p_penalties.gap_opening1 = gap_open1;   // O1 >= 0
	attributes.affine2p_penalties.gap_extension1 = gap_extension1; // E1 > 0
	attributes.affine2p_penalties.gap_opening2 = gap_open2;  // O2 >= 0
	attributes.affine2p_penalties.gap_extension2 = gap_extension2; // E2 > 0

	attributes.memory_mode = wavefront_memory_med;

	//attributes.heuristic.strategy = wf_heuristic_wfadaptive;
	//attributes.heuristic.min_wavefront_length = 10;
	//attributes.heuristic.max_distance_threshold = 50;
	//attributes.heuristic.steps_between_cutoffs = 1;

}

// Function to align two sequences based on given rare match pairs (anchors) and save the results.
void PairAligner::alignPairSeq(const std::vector<SequenceInfo>& data, RareMatchPairs anchors) {
	// Define the whole sequence interval for both sequences.
	Interval interval(0, data[0].seq_len, 0, data[1].seq_len);

	// Calculate the length of the first sequence.
	uint_t fst_length = data[0].seq_len;

	// Convert the rare match pairs (anchors) to intervals that need alignment.
	Intervals intervals_need_align = AnchorFinder::RareMatchPairs2Intervals(anchors, interval, fst_length);

	// Save the intervals that need alignment to a CSV file for further analysis or debugging.
	saveIntervalsToCSV(intervals_need_align, joinPaths(save_file_path, INTERVAL_NAME));

	// Perform the alignment on the intervals and return the resulting CIGAR string.
	cigar final_cigar = alignIntervals(data, intervals_need_align, anchors);

	// Save the resulting CIGAR string to a text file.
	saveCigarToTxt(final_cigar, joinPaths(save_file_path, CIGAR_NAME));

	// Optionally, verify the correctness of the generated CIGAR string against the input sequences.
	// verifyCigar(final_cigar, data);

	// Convert the final CIGAR string to a FASTA format and save it to a file for visualization or further analysis.
	cigarToFasta(final_cigar, data, joinPaths(save_file_path, FASTA_NAME));

	return; // End of the function.
}


// Function to align specified intervals within sequences and combine the result with anchor alignments.
cigar PairAligner::alignIntervals(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align, const RareMatchPairs& anchors) {
	// Initialize a vector to store aligned intervals as CIGAR strings for each interval.
	cigars aligned_interval_cigar(intervals_need_align.size());

	// Vector to keep track of which intervals actually need wavefront alignment.
	std::vector<uint_t> aligned_intervals_index;

	// Loop through each interval that needs alignment.
	for (uint_t i = 0; i < intervals_need_align.size(); ++i) {
		// Retrieve the current interval and its lengths for both sequences.
		Interval tmp_interval = intervals_need_align[i];
		uint_t fst_len = tmp_interval.len1;
		uint_t scd_len = tmp_interval.len2;
		// Extract the corresponding subsequences from both sequences.
		std::string seq1 = data[0].sequence.substr(tmp_interval.pos1, tmp_interval.len1);
		std::string seq2 = data[1].sequence.substr(tmp_interval.pos2, tmp_interval.len2);

		// Handle cases where one of the subsequences is empty.
		if (fst_len == 0) {
			aligned_interval_cigar[i] = cigar(1, cigarToInt('I', scd_len));
			continue;
		}
		if (scd_len == 0) {
			aligned_interval_cigar[i] = cigar(1, cigarToInt('D', fst_len));
			continue;
		}
		// Handle specific conditions based on the lengths of the subsequences.
		if (fst_len <= 5 && scd_len > 100) {
			cigar tmp_cigar;
			uint_t cigar_len = 1;
			char cur_state = (seq1[0] == seq2[0]) ? '=' : 'X';

			// Compare each character and construct the CIGAR string accordingly.
			for (uint_t i = 1; i < fst_len; ++i) {
				if (seq1[i] == seq2[i]) {
					if (cur_state == '=') {
						++cigar_len;
					}
					else {
						tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));
						cur_state = '=';
						cigar_len = 1;
					}
				}
				else {
					if (cur_state == 'X') {
						++cigar_len;
					}
					else {
						tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));
						cur_state = 'X';
						cigar_len = 1;
					}
				}
			}
			// Add the last operation to the CIGAR string.
			tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));

			// Handle the case where the second sequence is longer.
			if (scd_len > fst_len) {
				tmp_cigar.emplace_back(cigarToInt('I', scd_len - fst_len));
			}

			aligned_interval_cigar[i] = tmp_cigar;
			continue;
		}

		// Similar handling for the opposite condition.
		if (scd_len <= 5 && fst_len > 100) {
			cigar tmp_cigar;
			uint_t cigar_len = 1;
			char cur_state = (seq1[0] == seq2[0]) ? '=' : 'X';

			for (uint_t i = 1; i < scd_len; ++i) {
				if (seq1[i] == seq2[i]) {
					if (cur_state == '=') {
						++cigar_len;
					}
					else {
						tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));
						cur_state = '=';
						cigar_len = 1;
					}
				}
				else {
					if (cur_state == 'X') {
						++cigar_len;
					}
					else {
						tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));
						cur_state = 'X';
						cigar_len = 1;
					}
				}
			}

			tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));

			if (fst_len > scd_len) {
				tmp_cigar.emplace_back(cigarToInt('D', fst_len - scd_len));
			}

			aligned_interval_cigar[i] = tmp_cigar;
			continue;
		}

		// Add the index to the list for wavefront alignment if it wasn't handled above.
		aligned_intervals_index.emplace_back(i);
	}
	// Perform wavefront alignment for the intervals needing it.
	alignIntervalsUsingWavefront(data, intervals_need_align, aligned_intervals_index, aligned_interval_cigar);

	// Print debug information for the aligned intervals.
	printCigarDebug(data, aligned_interval_cigar, intervals_need_align);

	// Combine the aligned intervals with anchor alignments and return the final CIGAR string.
	return combineCigarsWithAnchors(aligned_interval_cigar, anchors);
}



void PairAligner::verifyCigar(const cigar& final_cigar, const std::vector<SequenceInfo>& data) {
	if (data.size() < 2) {
		logger.error() << "Not enough sequences provided for verification.\n";
		return;
	}

	const std::string& pattern = data[0].sequence;
	const std::string& text = data[1].sequence;
	uint_t pattern_pos = 0, text_pos = 0;

	for (const auto& unit : final_cigar) {
		char operation;
		uint32_t len;

		intToCigar(unit, operation, len); // Convert each cigar unit back to operation and length

		switch (operation) {
		case '=': // Sequence match
			for (uint_t i = 0; i < len; ++i) {
				if (pattern[pattern_pos + i] != text[text_pos + i]) {
					logger.error() << "Mismatch found where exact match expected at seq1 position "
						<< pattern_pos + i << " and seq2 position " << text_pos + i << std::endl;
					return;
				}
			}
			pattern_pos += len;
			text_pos += len;
			break;
		case 'X': // Mismatch
			for (uint_t i = 0; i < len; ++i) {
				if (pattern[pattern_pos + i] == text[text_pos + i]) {
					logger.error() << "Exact match found where mismatch expected at seq1 position "
						<< pattern_pos + i << " and seq2 position " << text_pos + i << std::endl;
					return;
				}
			}
			pattern_pos += len;
			text_pos += len;
			break;
		case 'M': // Generic match/mismatch
			pattern_pos += len;
			text_pos += len;
			break;
		case 'I': // Insertion
			text_pos += len;
			break;
		case 'D': // Deletion
			pattern_pos += len;
			break;
		default:
			std::cerr << "Unknown CIGAR operation '" << operation << "' encountered.\n";
			return;
		}
	}

	// Check if the end positions match the sequence lengths
	if (pattern_pos != data[0].seq_len || text_pos != data[1].seq_len) {
		logger.error() << "CIGAR does not fully align sequences. Seq1 aligned length: " << pattern_pos
			<< ", Seq2 aligned length: " << text_pos << std::endl;
	}
	else {
		logger.info() << "CIGAR verification successful.\n";
	}
}

// Saves the final CIGAR string to a text file.
void PairAligner::saveCigarToTxt(const cigar& final_cigar, const std::string& filename) {
	// Open the specified file for writing.
	std::ofstream outFile(filename);
	// Check if the file was successfully opened.
	if (!outFile.is_open()) {
		// Log an error message if the file could not be opened.
		logger.error() << "Error: Unable to open file " << filename << " for writing.\n";
		return; // Exit the function if the file cannot be opened.
	}

	// Iterate through each unit in the final CIGAR string.
	for (const auto& unit : final_cigar) {
		char operation; // Variable to hold the operation character.
		uint32_t len; // Variable to hold the length of the operation.
		// Convert the unit from its compact integer form back to operation and length.
		intToCigar(unit, operation, len);

		// Write the length and operation to the file.
		outFile << len << operation;
	}

	// End the CIGAR string with a newline character.
	outFile << std::endl;
	// Close the file after writing.
	outFile.close();

	// Log a message indicating the CIGAR string has been successfully saved.
	logger.info() << "CIGAR has been saved to " << filename << std::endl;
}


void PairAligner::printCigarDebug(const std::vector<SequenceInfo>& data, const cigars& aligned_interval_cigar, const Intervals& intervals_need_align) {
	for (uint_t i = 0; i < aligned_interval_cigar.size(); i++) {
		Interval tmp_interval = intervals_need_align[i];
		std::string seq1 = data[0].sequence.substr(tmp_interval.pos1, tmp_interval.len1);
		std::string seq2 = data[1].sequence.substr(tmp_interval.pos2, tmp_interval.len2);
		logger.debug() << "CIGAR: " << i + 1 << "\n";
		logger.debug() << "\n" << seq1 << "\n" << seq2 << "\n";
		for (uint_t j = 0; j < aligned_interval_cigar[i].size(); j++) {
			char operation;
			uint32_t len;
			intToCigar(aligned_interval_cigar[i][j], operation, len);
			logger.debug() << operation << len << "\n";
		}
	}
}

// This function combines multiple cigar vectors and intersperses them with 'anchor' operations. 
// It also removes any zero-length operations from the start and end of the combined vector.
cigar PairAligner::combineCigarsWithAnchors(const cigars& aligned_interval_cigar, const RareMatchPairs& anchors) {
	cigar final_cigar;
	auto anchor_cigar = anchors.begin();

	std::string confidence_csv = joinPaths(save_file_path, CONFIDENCE_CSV);

	std::ofstream csv_file(confidence_csv);
	if (!csv_file.is_open()) {
		logger.error() << "Error opening file " << confidence_csv << std::endl;
	}
	csv_file << "cigar,confidence,rare match\n";

	// Iterate through each cigar vector and add its units to the final_cigar vector.
	// Also intersperse 'anchor' operations between the cigar vectors.
	for (const auto& single_cigar : aligned_interval_cigar) {
		if (single_cigar.size() == 1) {
			char operation;
			uint32_t len;
			intToCigar(single_cigar[0], operation, len);
			if (len > 0) {
				csv_file << len << operation << "," << 1 << "," << 0 << "\n";
				final_cigar.push_back(single_cigar[0]);
			}

		}
		else {
			for (const auto& unit : single_cigar) {
				char operation;
				uint32_t len;
				intToCigar(unit, operation, len);
				if (len > 0) {
					csv_file << len << operation << "," << 0 << "," << 0 << "\n";
					final_cigar.push_back(unit);
				}
			}
		}

		// If there's an anchor, add an '=' operation with its match_length.
		if (anchor_cigar != anchors.end()) {
			final_cigar.push_back(cigarToInt('=', anchor_cigar->match_length));
			csv_file << anchor_cigar->match_length << "=," << 1 << "," << 1 << "\n";
			++anchor_cigar;
		}

	}

	// Remove zero-length operations from the start of the final_cigar, if present.
	if (!final_cigar.empty()) {
		char operation;
		uint32_t len;
		intToCigar(final_cigar.front(), operation, len);
		if (len == 0) {
			final_cigar.erase(final_cigar.begin());
		}
	}

	// Remove zero-length operations from the end of the final_cigar, if present.
	if (!final_cigar.empty()) {
		char operation;
		uint32_t len;
		intToCigar(final_cigar.back(), operation, len);
		if (len == 0) {
			final_cigar.pop_back();
		}
	}
	csv_file.close();
	return final_cigar;
}

void PairAligner::cigarToFasta(const cigar& final_cigar, const std::vector<SequenceInfo>& data, const std::string& fasta_filename) {
	if (data.size() < 2) {
		logger.error() << "Not enough sequences provided for CIGAR to FASTA conversion.\n";
		return;
	}

	// Open the output FASTA file.
	std::ofstream fasta_file(fasta_filename);
	if (!fasta_file.is_open()) {
		logger.error() << "Failed to open FASTA output file: " << fasta_filename << std::endl;
		return;
	}

	const std::string& pattern = data[0].sequence;
	const std::string& text = data[1].sequence;
	uint_t pattern_pos = 0, text_pos = 0;
	std::string aligned_seq1, aligned_seq2;

	// Process each unit in the CIGAR string to construct the aligned sequences.
	for (const auto& unit : final_cigar) {
		char operation;
		uint32_t len;

		intToCigar(unit, operation, len); // Convert each cigar unit back to operation and length.

		switch (operation) {
		case '=': // Sequence match
		case 'X': // Mismatch
		case 'M': // Generic match/mismatch
			aligned_seq1.append(pattern.substr(pattern_pos, len));
			aligned_seq2.append(text.substr(text_pos, len));
			pattern_pos += len;
			text_pos += len;
			break;
		case 'I': // Insertion
			aligned_seq1.append(len, '-'); // Add gaps to sequence1
			aligned_seq2.append(text.substr(text_pos, len));
			text_pos += len;
			break;
		case 'D': // Deletion
			aligned_seq1.append(pattern.substr(pattern_pos, len));
			aligned_seq2.append(len, '-'); // Add gaps to sequence2
			pattern_pos += len;
			break;
		default:
			logger.error() << "Unknown CIGAR operation '" << operation << "' encountered.\n";
			fasta_file.close();
			return;
		}
	}

	fasta_file << ">" << data[0].header << std::endl;
	fasta_file << aligned_seq1 << "\n";

	fasta_file << ">" << data[1].header << std::endl;
	fasta_file << aligned_seq2 << "\n";

	fasta_file.close();
	logger.info() << fasta_filename << " has been saved successfully!" << std::endl;
}

// Function to align sequences within specified intervals using the wavefront alignment method.
void PairAligner::alignIntervalsUsingWavefront(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align, std::vector<uint_t>& aligned_intervals_index, cigars& aligned_interval_cigar) {
	ThreadPool pool(thread_num); // Create a thread pool with the determined number of threads.
	logger.info() << "Begin to align intervals using wavefront alignment method." << std::endl;

	// Iterate through each interval that requires alignment.
	for (uint_t i = 0; i < aligned_intervals_index.size(); ++i) {
		uint_t index = aligned_intervals_index[i]; // Get the index of the current interval.
		Interval tmp_interval = intervals_need_align[index]; // Retrieve the interval details.
		// Extract the subsequences from both sequences based on the interval information.
		std::string seq1 = data[0].sequence.substr(tmp_interval.pos1, tmp_interval.len1);
		std::string seq2 = data[1].sequence.substr(tmp_interval.pos2, tmp_interval.len2);

		// Check if parallel processing is enabled.
		if (thread_num) {
			// If parallel processing is enabled, enqueue alignment tasks to the thread pool.
			pool.enqueue([this, seq1, seq2, &aligned_interval_cigar, index]() {
				// Create a new wavefront aligner instance with specified attributes.
				wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
				// Perform the alignment using the wavefront aligner.
				wavefront_align(wf_aligner, seq1.c_str(), seq1.length(), seq2.c_str(), seq2.length());
				uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
				int cigar_length; // Length of the CIGAR string.
				// Retrieve the CIGAR string from the wavefront aligner.
				cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
				// Convert the CIGAR buffer to a vector and store it in the aligned_interval_cigar vector.
				aligned_interval_cigar[index] = convertToCigarVector(cigar_buffer, cigar_length);
				wavefront_aligner_delete(wf_aligner); // Free the aligner resources. 
				});
		}
		else {
			// If parallel processing is not enabled, perform the alignment in the main thread.
			wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
			wavefront_align(wf_aligner, seq1.c_str(), seq1.length(), seq2.c_str(), seq2.length());
			uint32_t* cigar_buffer;
			int cigar_length;
			cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
			aligned_interval_cigar[index] = convertToCigarVector(cigar_buffer, cigar_length);
			wavefront_aligner_delete(wf_aligner); // Free the aligner resources.
		}
	}

	if (thread_num) {
		pool.waitAllTasksDone(); // Wait for all alignment tasks in the thread pool to complete.
	}
	logger.info() << "Wavefront alignment of intervals has been completed." << std::endl;
}





