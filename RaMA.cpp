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
#include "gsacak.h"
#include "utils.h"
#include "logging.h"
#include "anchor.h"
#include "threadpool.h"
#include "pairwise_alignment.h"
#include "argparser.h"
extern "C" {
#include "wavefront/wavefront_align.h"
}
#include "Alignment/WFA2-lib/bindings/cpp/WFAligner.hpp"

Logger logger("/mnt/f/code/vs_code/RaMA/output/", "RaMA", true, info);


int main(int argc, char** argv) {
    logger.info() << "Start RaMA!" << std::endl;
    std::ios::sync_with_stdio(false);


    Parser p(argc, argv);
    p.add("-i", "--input", "The path of input fasta data.", Mode::REQUIRED);
    p.add("-o", "--output", "The path of output directory.", Mode::REQUIRED);

    p.add("-t", "--threads", "The number of threads to use.", Mode::OPTIONAL);

    p.add("-m", "--match", "The match score for sequence alignment.", Mode::OPTIONAL);
    p.add("-s", "--mismatch", "The mismatch score for sequence alignment.", Mode::OPTIONAL);
    p.add("-g", "--gap_open1", "The gap open score for sequence alignment.", Mode::OPTIONAL);
    p.add("-e", "--gap_extension1", "The gap extension score for sequence alignment.", Mode::OPTIONAL);
    p.add("-G", "--gap_open2", "The gap open score for sequence alignment.", Mode::OPTIONAL);
    p.add("-E", "--gap_extension2", "The gap extension score for sequence alignment.", Mode::OPTIONAL);


    auto args = p.parse();

    if (!args.parsedSuccessfully()) {
        logger.error() << "Unsuccessful parse args\n";
        p.printHelpString();
        logger.info() << "Exit RaMA!" << std::endl;
        return -1;
    }

	
	const char* data_path = "/mnt/f/code/vs_code/RaMA/data/human.fasta";
    RareMatchPairs final_anchors;
	std::vector<SequenceInfo>* data = new std::vector<SequenceInfo>(readDataPath(data_path));
    {
        AnchorFinder anchor_finder(*data, 10, "/mnt/f/code/vs_code/RaMA/output/save/", false, false);
        final_anchors = anchor_finder.lanuchAnchorSearching();
    }

    PairAligner pair_aligner(0, 3, 4, 2, 12, 1, 10);
    pair_aligner.alignPairSeq(*data, final_anchors);
    logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;
	delete data;
    logger.info() << "End RaMA!" << std::endl;

    return 0;

}