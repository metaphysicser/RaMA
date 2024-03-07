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

PairAligner::PairAligner(int_t match, int_t mismatch, int_t gap_open, int_t gap_extension, bool use_parallel): 
	match(match), 
	mismatch(mismatch),
	gap_open(gap_open),
	gap_extension(gap_extension),
	use_parallel(use_parallel) {
	score_affine = seqan2::Score<int16_t, seqan2::Simple>(match, mismatch, gap_open, gap_extension);
	setNumThreads(exec_policy, std::thread::hardware_concurrency());
}

void PairAligner::alignPairSeq(const std::vector<SequenceInfo>& data, RareMatchPairs anchors) {
	Interval first_interval = { 0, data[0].seq_len };
	Interval second_interval = { 0, data[1].seq_len };
	uint_t fst_length = data[0].seq_len;
	Intervals intervals_need_align = AnchorFinder::RareMatchPairs2Intervals(anchors, first_interval, second_interval, fst_length);
	saveIntervalsToCSV(intervals_need_align, "/mnt/f/code/vs_code/RaMA/output/intervals_need_align.csv");
	alignIntervals(data, intervals_need_align);
	return;
}


void PairAligner::alignIntervals(const std::vector<SequenceInfo>& data, const Intervals& intervals_need_align) {
	std::vector<std::vector<std::string>> aligned_intervals(2);
	for (auto& seq : aligned_intervals) {
		seq.resize(intervals_need_align.size());
	}
	std::vector<uint_t> aligned_intervals_index;

	seqan2::StringSet<TAlignedSequence> seqs1;
	seqan2::StringSet<TAlignedSequence> seqs2;


	for (uint_t i = 0; i < intervals_need_align.size(); ++i) {
		std::pair<Interval, Interval> tmp_interval_pairs = intervals_need_align[i];
		if (tmp_interval_pairs.first.second == 0) {
			aligned_intervals[0][i] = std::string(tmp_interval_pairs.second.second, '-');
			aligned_intervals[1][i] = data[1].sequence.substr(tmp_interval_pairs.second.first, tmp_interval_pairs.second.second);
			continue;
		}

		if (tmp_interval_pairs.second.second == 0) {
			aligned_intervals[0][i] = data[0].sequence.substr(tmp_interval_pairs.first.first, tmp_interval_pairs.first.second);
			aligned_intervals[1][i] = std::string(tmp_interval_pairs.first.second, '-');
			continue;
		}

		TSequence seq1 = data[0].sequence.substr(tmp_interval_pairs.first.first, tmp_interval_pairs.first.second).c_str();
		TSequence seq2 = data[1].sequence.substr(tmp_interval_pairs.second.first, tmp_interval_pairs.second.second).c_str();

		appendValue(seqs1, TAlignedSequence(seq1));
		appendValue(seqs2, TAlignedSequence(seq2));

		aligned_intervals_index.push_back(i);
	
	}

	seqan2::String<int16_t> scores = seqan2::globalAlignment(exec_policy, seqs1, seqs2, score_affine);
	std::cout << seqan2::length(seqs1) << std::endl;
	//for (int16_t score : scores)
	//    std::cout << "Score: " << score << "\n";

	for (size_t pos = 0; pos < seqan2::length(seqs1); ++pos) // print out alignments
	{
		
		std::cout << seqs1[pos] << "\n";
		std::cout << seqs2[pos] << "\n\n";
	}
	return;
}

