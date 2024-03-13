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

PairAligner::PairAligner(int_t match, int_t mismatch, int_t gap_open1, int_t gap_extension1, int_t gap_open2, int_t gap_extension2, bool use_parallel):
	match(match), 
	mismatch(mismatch),
	gap_open1(gap_open1),
	gap_extension1(gap_extension1),
	gap_open2(gap_open2),
	gap_extension2(gap_extension2),
	use_parallel(use_parallel) {
	attributes = wavefront_aligner_attr_default;
	attributes.distance_metric = gap_affine_2p;
	attributes.affine2p_penalties.match = match;
	attributes.affine2p_penalties.mismatch = mismatch;       // X > 0
	attributes.affine2p_penalties.gap_opening1 = gap_open1;   // O1 >= 0
	attributes.affine2p_penalties.gap_extension1 = gap_extension1; // E1 > 0
	attributes.affine2p_penalties.gap_opening2 = gap_open2;  // O2 >= 0
	attributes.affine2p_penalties.gap_extension2 = gap_extension2; // E2 > 0
	/*attributes.alignment_scope = compute_alignment;
	attributes.alignment_form.span = alignment_end2end;
	attributes.heuristic.strategy = wf_heuristic_wfadaptive;
	attributes.heuristic.min_wavefront_length = 10;
	attributes.heuristic.max_distance_threshold = 50;
	attributes.heuristic.steps_between_cutoffs = 1;*/
}

void PairAligner::alignPairSeq(const std::vector<SequenceInfo>& data, RareMatchPairs anchors) {
	Interval first_interval = { 0, data[0].seq_len };
	Interval second_interval = { 0, data[1].seq_len };
	uint_t fst_length = data[0].seq_len;
	Intervals intervals_need_align = AnchorFinder::RareMatchPairs2Intervals(anchors, first_interval, second_interval, fst_length);
	saveIntervalsToCSV(intervals_need_align, "/mnt/f/code/vs_code/RaMA/output/intervals_need_align.csv");
	cigar final_cigar = alignIntervals(data, intervals_need_align, anchors);
	saveCigarToTxt(final_cigar, "/mnt/f/code/vs_code/RaMA/output/final_cigar.txt");
	verifyCigar(final_cigar, data);
	return;
}


cigar PairAligner::alignIntervals(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align, const RareMatchPairs& anchors) {
	cigars aligned_interval_cigar(intervals_need_align.size());
	
	std::vector<uint_t> aligned_intervals_index;

	// todo: tmp_interval_pairs.first.second change name
	for (uint_t i = 0; i < intervals_need_align.size(); ++i) {
		std::pair<Interval, Interval> tmp_interval_pairs = intervals_need_align[i];
		uint_t fst_len = tmp_interval_pairs.first.second;
		uint_t scd_len = tmp_interval_pairs.second.second;
		std::string seq1 = data[0].sequence.substr(tmp_interval_pairs.first.first, tmp_interval_pairs.first.second);
		std::string seq2 = data[1].sequence.substr(tmp_interval_pairs.second.first, tmp_interval_pairs.second.second);

		if (fst_len == 0) { 
			// Insert
			aligned_interval_cigar[i] = cigar(1, cigarToInt('I', scd_len));
			continue;
		}

		if (scd_len == 0) {
			aligned_interval_cigar[i] = cigar(1, cigarToInt('D', fst_len));
			continue;
		}
		// todo：这里的参数设置可以根据罚分调整
		if (fst_len <= 5 && scd_len > 100) {
			cigar tmp_cigar;
			uint_t cigar_len = 1; 
			char cur_state = (seq1[0] == seq2[0]) ? '=' : 'X'; 

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

			tmp_cigar.emplace_back(cigarToInt(cur_state, cigar_len));

			if (scd_len > fst_len) {
				tmp_cigar.emplace_back(cigarToInt('I', scd_len - fst_len));
			}

			aligned_interval_cigar[i] = tmp_cigar;
			continue;
		}


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

		

		aligned_intervals_index.emplace_back(i);

	}
	uint_t num_threads = std::thread::hardware_concurrency();
	ThreadPool pool(num_threads);

	for (uint_t i = 0; i < aligned_intervals_index.size(); ++i) {
		uint_t index = aligned_intervals_index[i];
		std::pair<Interval, Interval> tmp_interval_pairs = intervals_need_align[index];
		std::string seq1 = data[0].sequence.substr(tmp_interval_pairs.first.first, tmp_interval_pairs.first.second);
		std::string seq2 = data[1].sequence.substr(tmp_interval_pairs.second.first, tmp_interval_pairs.second.second);
		if (use_parallel) {
			pool.enqueue([this, seq1, seq2, &aligned_interval_cigar, index]() {
			wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
			const char* pattern = seq1.c_str();
			const char* text = seq2.c_str();
			wavefront_align(wf_aligner, pattern, strlen(pattern), text, strlen(text)); // Align
			uint32_t* cigar_buffer;
			int cigar_length;
			cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
			aligned_interval_cigar[index] = convertToCigarVector(cigar_buffer, cigar_length);
			wavefront_aligner_delete(wf_aligner); // Free	
			});
		} else {
			wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
			const char* pattern = seq1.c_str();
			const char* text = seq2.c_str();
			wavefront_align(wf_aligner, pattern, strlen(pattern), text, strlen(text)); // Align
			uint32_t* cigar_buffer;
			int cigar_length;
			cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
			aligned_interval_cigar[index] = convertToCigarVector(cigar_buffer, cigar_length);
			wavefront_aligner_delete(wf_aligner); // Free			
		}
		
	}	
	if(use_parallel)
		pool.waitAllTasksDone(); // Wait for all tasks to complete

	for (uint_t i = 0; i < aligned_interval_cigar.size(); i++) {
		std::pair<Interval, Interval> tmp_interval_pairs = intervals_need_align[i];
		std::string seq1 = data[0].sequence.substr(tmp_interval_pairs.first.first, tmp_interval_pairs.first.second);
		std::string seq2 = data[1].sequence.substr(tmp_interval_pairs.second.first, tmp_interval_pairs.second.second);
		logger.debug() << "CIGAR: " << i+1  << "\n";
		logger.debug() << "\n" << seq1 << "\n" << seq2 << "\n";
		for (uint_t j = 0; j < aligned_interval_cigar[i].size(); j++) {
			char operation;
			uint32_t len;
			intToCigar(aligned_interval_cigar[i][j], operation, len);
			logger.debug() << operation << len << "\n";
		}
	}

	cigar final_cigar;
	auto anchor_cigar = anchors.begin();
	for (const auto& single_cigar : aligned_interval_cigar) {
		for (const auto& unit : single_cigar) {
			final_cigar.push_back(unit);
		}
		if (anchor_cigar != anchors.end()) {
			final_cigar.push_back(cigarToInt('=', anchor_cigar->match_length));
			anchor_cigar++;
		}
	}

	char operation;
	uint_t len;

	intToCigar(final_cigar.front(), operation, len);
	if (len == 0) {
		final_cigar.erase(final_cigar.begin());
	}

	if (!final_cigar.empty()) {
		intToCigar(final_cigar.back(), operation, len);
		if (len == 0) {
			final_cigar.pop_back();
		}
	}
	return final_cigar;
}


cigar PairAligner::convertToCigarVector(uint32_t* cigar_buffer, int cigar_length) {
	cigar resultCigar;
	for (int i = 0; i < cigar_length; ++i) {
		resultCigar.push_back(cigar_buffer[i]);
	}
	return resultCigar;
}

uint32_t PairAligner::cigarToInt(char operation, uint32_t len) {
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


void PairAligner::intToCigar(uint32_t cigar, char& operation, uint32_t& len) {
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

void PairAligner::verifyCigar(const cigar& final_cigar, const std::vector<SequenceInfo>& data) {
	if (data.size() < 2) {
		logger.error() << "Error: Not enough sequences provided for verification.\n";
		return;
	}

	const std::string& pattern = data[0].sequence;
	const std::string& text = data[1].sequence;
	uint_t pattern_pos = 0, text_pos = 0;

	for (const auto& unit : final_cigar) {
		char operation;
		uint_t len;
		intToCigar(unit, operation, len); // Convert each cigar unit back to operation and length
		switch (operation) {
		case '=': // Sequence match
			for (uint_t i = 0; i < len; ++i) {
				if (pattern[pattern_pos + i] != text[text_pos + i]) {
					logger.error() << "Mismatch found where exact match expected at pattern position "
						<< pattern_pos + i << " and text position " << text_pos + i << std::endl;
					return;
				}
			}
			pattern_pos += len;
			text_pos += len;
			break;
		case 'X': // Mismatch
			for (uint_t i = 0; i < len; ++i) {
				if (pattern[pattern_pos + i] == text[text_pos + i]) {
					std::cerr << "Exact match found where mismatch expected at pattern position "
						<< pattern_pos + i << " and text position " << text_pos + i << std::endl;
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
		logger.error() << "CIGAR does not fully align sequences. Pattern aligned length: " << pattern_pos
			<< ", Text aligned length: " << text_pos << std::endl;
	}
	else {
		logger.info() << "CIGAR verification successful.\n";
	}
}

void PairAligner::saveCigarToTxt(const cigar& final_cigar, const std::string& filename) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		logger.error() << "Error: Unable to open file " << filename << " for writing.\n";
		return;
	}

	for (const auto& unit : final_cigar) {
		char operation;
		uint32_t len;
		intToCigar(unit, operation, len);

		outFile << len << operation;
	}

	outFile << std::endl;
	outFile.close();

	logger.info() << "CIGAR has been saved to " << filename << std::endl;
}



