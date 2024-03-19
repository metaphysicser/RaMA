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

Logger logger("RaMA", true, info);

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);

    Parser p(argc, argv);
    p.add("-i", "--input", "Specify the path to the input FASTA file. This file should contain the DNA you wish to align.", Mode::REQUIRED);
    p.add("-o", "--output", "Define the path to the output directory where all generated files will be saved. This includes alignment results, anchor binary files, and any other output.", Mode::REQUIRED);

    p.add("-t", "--threads", "Set the number of threads to use during the alignment process. If not specified, the program will use 0 as a default value.", Mode::OPTIONAL);

    p.add("-s", "--save", "Enable this flag to save the anchor binary file to the specified output directory. This can be useful for reusing the SA, LCP, Linear Sparse Table and so on in future alignments without recalculating them.", Mode::BOOLEAN);
    p.add("-l", "--load", "Enable this flag to load an existing anchor binary file from the specified output directory. This skips the anchor SA, LCP and Linear Sparse Table construction and directly uses the them for finding anchors.", Mode::BOOLEAN);

    p.add("-c", "--max_match_count", "Define the maximum number of rare matches to consider when finding anchors. We don't recommend you to reset this value.", Mode::OPTIONAL);

    p.add("-m", "--match", "Sets the score for matching bases in sequence alignment. A higher value increases the incentive for aligning matching characters, enhancing alignment accuracy. This parameter plays a crucial role in balancing between matches and gaps, especially in the context of dual-cost gap-affine distances where the treatment of gaps varies based on their length. Default value is 0.", Mode::OPTIONAL);
    p.add("-x", "--mismatch", "Defines the penalty for mismatching bases. In the framework of dual-cost gap-affine distances, a higher mismatch penalty disincentivizes the alignment from incorporating mismatches, which is particularly important when considering the alignment's overall strategy towards managing short and long gaps. Default value is 3.", Mode::OPTIONAL);
    p.add("-g", "--gap_open1", "Specifies the penalty for initiating a short gap in the alignment. This penalty is a key component of the dual-cost gap-affine distance, allowing for a differentiated approach to short and long gaps. A well-calibrated short gap opening penalty can prevent unnecessary fragmentation of the alignment by sporadic mismatches. Default value is 4, reflecting its role as a penalty.", Mode::OPTIONAL);
    p.add("-e", "--gap_extension1", "Determines the penalty for extending an already existing short gap. This parameter complements the short gap opening penalty by influencing the cost of continuing a gap once it has been opened. Adjusting this penalty allows for fine-tuning the alignment's sensitivity to continuous versus sporadic gaps. Default value is 2, facilitating the extension of gaps that have already been penalized for initiating.", Mode::OPTIONAL);
    p.add("-G", "--gap_open2", "Similar to --gap_open1 but applies to the initiation of long gaps. In the context of dual-cost gap-affine distances, this penalty differentiates between the costs associated with starting short versus long gaps, acknowledging the distinct challenges posed by longer gaps in sequence alignment. This parameter allows for the strategic management of long gaps, aiming to reduce the fragmentation of the alignment. Default value is 12.", Mode::OPTIONAL);
    p.add("-E", "--gap_extension2", "Mirrors --gap_extension1 but for extending long gaps. This parameter is crucial for managing the continuity of long gaps once they have been opened, playing a significant role in the dual-cost gap-affine distance approach. By setting a specific penalty for long gap extensions, it provides a mechanism for handling long, continuous gaps more leniently than short gaps, reflecting their different impacts on alignment quality. Default value is 1.", Mode::OPTIONAL);

    auto args = p.parse();

    if (!args.parsedSuccessfully()) {
        std::cerr << "Unsuccessful parse args\n";
        p.printHelpString();
        std::cout << "Exit RaMA!" << std::endl;
        return -1;
    }

    std::string data_path, output_path;
    bool save, load;
    uint_t thread_num, max_match_count;
    int_t match, mismatch, gap_open1, gap_open2, gap_extension1, gap_extension2;

    try {
        data_path = args["--input"];
        output_path = args["--output"];
        thread_num = args["--threads"].empty()? 0 : std::stoi(args["--threads"]);
        save = args["--save"] == "1";
        load = args["--load"] == "1";
        max_match_count = getMaxValue(args["--max_match_count"].empty()? 100 : std::stoi(args["--max_match_count"]), 2);
        match = args["--match"].empty()? 0 : std::stoi(args["--match"]);
        mismatch = args["--mismatch"].empty()? 3 : std::stoi(args["--mismatch"]);
        gap_open1 = args["--gap_open1"].empty()? 4 : std::stoi(args["--gap_open1"]);
        gap_extension1 = args["--gap_extenstion1"].empty()? 2 : std::stoi(args["--gap_extension1"]);
        gap_open2 = args["--gap_open2"].empty() ? 12 : std::stoi(args["--gap_open2"]);
        gap_extension2 = args["--gap_extension2"].empty()? 1 : std::stoi(args["--gap_extension2"]);
    }
    catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Unsuccessful parse args\n";
        p.printHelpString();
		std::cout << "Exit RaMA!" << std::endl;
		return -1;
	}

    logger.setDir(output_path);
    logger.info() << "Start RaMA!" << std::endl;
       
	// const char* data_path = "/mnt/f/code/vs_code/RaMA/data/human.fasta";
    RareMatchPairs final_anchors;
	std::vector<SequenceInfo>* data = new std::vector<SequenceInfo>(readDataPath(data_path.c_str()));
    {
        AnchorFinder anchor_finder(*data, thread_num, output_path.c_str(), load, save, max_match_count);
        final_anchors = anchor_finder.lanuchAnchorSearching();
    }

    PairAligner pair_aligner(output_path, match, mismatch, gap_open1, gap_extension1, gap_open2, gap_extension2, thread_num);
    pair_aligner.alignPairSeq(*data, final_anchors);
    logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;
	delete data;
    logger.info() << "End RaMA!" << std::endl;

    return 0;
}