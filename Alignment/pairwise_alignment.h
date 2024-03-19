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
#pragma once

#include "gsacak.h"
#include "logging.h"
#include "utils.h"
#include "anchor.h"
extern "C" {
#include "wavefront/wavefront_align.h"
}
#include "Alignment/WFA2-lib/bindings/cpp/WFAligner.hpp"

// Define types for handling CIGAR strings.
using cigarunit = uint32_t; // Represents a single operation in a CIGAR string.
using cigar = std::vector<cigarunit>; // Represents a CIGAR string.
using cigars = std::vector<cigar>; // Represents a collection of CIGAR strings.

// Convert a CIGAR operation and its length to a compact integer representation.
uint32_t cigarToInt(char operation, uint32_t len);

// Convert a compact integer representation of a CIGAR operation back to its character and length.
void intToCigar(uint32_t cigar, char& operation, uint32_t& len);

// Convert a buffer of compact integer CIGAR operations to a vector representation.
cigar convertToCigarVector(uint32_t* cigar_buffer, int cigar_length);

// Class for performing pairwise sequence alignment.
class PairAligner {
private:
	// Scoring parameters for sequence alignment.
	int_t match;
	int_t mismatch;
	int_t gap_open1;
	int_t gap_extension1;
	int_t gap_open2;
	int_t gap_extension2;

	uint_t thread_num; // Flag to enable parallel processing.

	wavefront_aligner_attr_t attributes; // Attributes for the wavefront aligner.

	// Align intervals within sequences and return the resulting CIGAR string.
	cigar alignIntervals(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align, const RareMatchPairs& anchors);

	// Verify the correctness of the generated CIGAR string against the input sequences.
	void verifyCigar(const cigar& final_cigar, const std::vector<SequenceInfo>& data);

	// Save a CIGAR string to a text file.
	void saveCigarToTxt(const cigar& final_cigar, const std::string& filename);

	// Print debug information for aligned interval CIGAR strings.
	void printCigarDebug(const std::vector<SequenceInfo>& data, const cigars& aligned_interval_cigar, const Intervals& intervals_need_align);

	// Combine individual CIGAR strings with anchor alignments.
	cigar combineCigarsWithAnchors(const cigars& aligned_interval_cigar, const RareMatchPairs& anchors);

	// Convert Cigar to fasta file.
	void cigarToFasta(const cigar& final_cigar, const std::vector<SequenceInfo>& data, const std::string& fasta_filename);

	// Use the wavefront alignment algorithm to align sequence intervals.
	void alignIntervalsUsingWavefront(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align, std::vector<uint_t>& aligned_intervals_index, cigars& aligned_interval_cigar);

public:
	// Constructor to initialize the PairAligner with scoring parameters and parallel processing flag.
	explicit PairAligner(int_t match = 0, int_t mismatch = 3, int_t gap_open1 = 4, int_t gap_extension1 = 2, int_t gap_open2 = 12, int_t gap_extension2 = 1, uint_t thread_num = 0);

	// Perform pairwise sequence alignment using provided data and optional anchors.
	void alignPairSeq(const std::vector<SequenceInfo>& data, RareMatchPairs anchors = {});

	// Overload of alignPairSeq to allow calling without explicitly specifying anchors.
	void alignPairSeq(const std::vector<SequenceInfo>& data) {
		alignPairSeq(data, {}); // Calling the first method with default second argument
	}
};
