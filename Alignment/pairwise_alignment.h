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
#include "config.h"
#include "logging.h"
#include "utils.h"
#include "anchor.h"
extern "C" {
#include "wavefront/wavefront_align.h"
}
#include "Alignment/WFA2-lib/bindings/cpp/WFAligner.hpp"

using cigarunit = uint32_t;
using cigar = std::vector<cigarunit>;
using cigars = std::vector<cigar>;

class PairAligner {
private:
	int_t match;
	int_t mismatch;
	int_t gap_open1;
	int_t gap_extension1;
	int_t gap_open2;
	int_t gap_extension2;

	bool use_parallel;

	wavefront_aligner_attr_t attributes;

	cigar alignIntervals(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align, const RareMatchPairs& anchors);

	void cigarToAlignment(cigar final_ciagr, const std::vector<SequenceInfo>& data);

	void verifyCigar(const cigar& final_cigar, const std::vector<SequenceInfo>& data);

	void saveCigarToTxt(const cigar& final_cigar, const std::string& filename);


public:
	explicit PairAligner(int_t match = 0, int_t mismatch = 3, int_t gap_open1 = 4, int_t gap_extension1 = 2, int_t gap_open2 = 12, int_t gap_extension2 = 1, bool use_parallel = true);

	void alignPairSeq(const std::vector<SequenceInfo>& data, RareMatchPairs anchors = {});

	void alignPairSeq(const std::vector<SequenceInfo>& data) {
		alignPairSeq(data, {}); // Calling the first method with default second argument
	}

	static uint32_t cigarToInt(char operation, uint32_t len);

	static void intToCigar(uint32_t cigar, char& operation, uint32_t& len);

	static cigar convertToCigarVector(uint32_t* cigar_buffer, int cigar_length);
};