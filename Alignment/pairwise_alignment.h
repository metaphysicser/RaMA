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
#include <seqan/align_parallel.h>


using TSequence = seqan2::String<seqan2::Dna>;
using TAlignedSequence = seqan2::Gaps<TSequence>;
using TThreadModel = seqan2::Parallel;
using TVectorSpec = seqan2::Vectorial;
using TExecPolicy = seqan2::ExecutionPolicy<TThreadModel, TVectorSpec>;

class PairAligner {
private:
	int_t match;
	int_t mismatch;
	int_t gap_open;
	int_t gap_extension;
	
	bool use_parallel;

	seqan2::Score<int16_t, seqan2::Simple> score_affine;
	TExecPolicy exec_policy;

	void alignIntervals(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align);

public:
	explicit PairAligner(int_t match, int_t mismatch, int_t gap_open, int_t gap_extension, bool use_parallel);

	void alignPairSeq(const std::vector<SequenceInfo>& data, RareMatchPairs anchors = {});

	void alignPairSeq(const std::vector<SequenceInfo>& data) {
		alignPairSeq(data, {}); // Calling the first method with default second argument
	}
};
