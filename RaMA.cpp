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
 // Created: 2024-01-29
#include "config.h"
#include "utils.h"
#include "logging.h"
#include "anchor.h"
#include "threadpool.h"
#include "pairwise_alignment.h"

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printint strings

#include <memory>
#include "D:/visualStudio/IDE/VC/Tools/MSVC/14.38.33130/include/memory"

Logger logger("/mnt/f/code/vs_code/RaMA/output/", "RaMA", true, info);

// Function to generate a random DNA sequence of given length
std::string generateRandomDNASequence(int length) {
    std::string bases = "ATCG";
    std::string sequence = "";
    for (int i = 0; i < length; ++i) {
        sequence += bases[rand() % 4]; // Randomly pick a base and append to the sequence
    }
    return sequence;
}

void f() {
    using TSequence = seqan2::String<seqan2::Dna>;
    using TAlignedSequence = seqan2::Gaps<TSequence>;
    using TThreadModel = seqan2::Parallel;
    using TVectorSpec = seqan2::Vectorial;
    using TExecPolicy = seqan2::ExecutionPolicy<TThreadModel, TVectorSpec>;

    TSequence seqH;
    TSequence seqV = generateRandomDNASequence(1000);

    seqan2::StringSet<TAlignedSequence> seqs1;
    seqan2::StringSet<TAlignedSequence> seqs2;

    /*resize(seqs1, 3);
    resize(seqs2, 3);

    assignValue(seqs1, 0, seqH);
    assignValue(seqs2, 0, seqH);

    assignValue(seqs1, 1, seqH);
    assignValue(seqs2, 1, seqH);

    assignValue(seqs1, 2, seqH);
    assignValue(seqs2, 2, seqH);*/

    //for (size_t i = 0; i < 1; ++i) // Create a data set of 100 dummy sequences
    //{
    //    appendValue(seqs1, TAlignedSequence(seqH));
    //    appendValue(seqs2, TAlignedSequence(seqV));
    //}



    //for (size_t i = 0; i < 2; ++i) // Create a data set of 100 dummy sequences
    //{
    //    appendValue(seqs1, TAlignedSequence(*new TSequence() = generateRandomDNASequence(10)));
    //    appendValue(seqs2, TAlignedSequence(*new TSequence() = generateRandomDNASequence(10)));
    //}

    for (size_t i = 0; i < 1; ++i) // Create a data set of 100 dummy sequences
    {
        appendValue(seqs1, TAlignedSequence(*new TSequence() = generateRandomDNASequence(20000000)));

        appendValue(seqs2, TAlignedSequence(*new TSequence() = generateRandomDNASequence(200000000)));

    }

    // arrayDestruct(seqs1.strings.data_begin, seqs1.strings.data_end);
    logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;


    /*TExecPolicy execPolicy;
    setNumThreads(execPolicy, 10);

    seqan2::Score<int16_t, seqan2::Simple> scoreAffine(0, -3, -2, -1);

    seqan2::String<int16_t> scores = seqan2::globalAlignment(execPolicy, seqs1, seqs2, scoreAffine);*/

   /* for (int16_t score : scores)
        std::cout << "Score: " << score << "\n";*/

    //for (size_t pos = 0; pos < seqan2::length(seqs1); ++pos) // print out alignments
    //{
    //    std::cout << seqs1[pos] << "\n";
    //    std::cout << seqs2[pos] << "\n\n";
    //}
}

int main(int argc, char** argv) {
	//std::ios::sync_with_stdio(false);
	//
	//const char* data_path = "/mnt/f/code/vs_code/RaMA/data/human.fasta";

	//std::vector<SequenceInfo>* data = new std::vector<SequenceInfo>(readDataPath(data_path));
	//AnchorFinder anchorfinder(*data, true,"/mnt/f/code/vs_code/RaMA/output/save/", false, true);
	//RareMatchPairs final_anchors = anchorfinder.lanuchAnchorSearching();
 //   PairAligner pair_aligner(0, -2, -3, -1, true);
 //   pair_aligner.alignPairSeq(*data, final_anchors);

	//delete data;
	logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;

    f();

    logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;





    return 0;
}