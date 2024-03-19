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
    Parser p(argc, argv);
    p.add("-i", "--input", "The path of input fasta data.", Mode::REQUIRED);
    p.add("-o", "--output", "The path of output directory.", Mode::REQUIRED);

    p.add("-t", "--threads", "The number of threads to use.", Mode::OPTIONAL);

    p.add("-s", "--save", "Save the Anchor binary file to ouput directory.", Mode::BOOLEAN);
    p.add("-l", "--load", "Load the Anchor binary file from ouput directory.", Mode::BOOLEAN);

    p.add("-m", "--match", "The match score for sequence alignment.", Mode::OPTIONAL);
    p.add("-x", "--mismatch", "The mismatch score for sequence alignment.", Mode::OPTIONAL);
    p.add("-g", "--gap_open1", "The gap open score for sequence alignment.", Mode::OPTIONAL);
    p.add("-e", "--gap_extension1", "The gap extension score for sequence alignment.", Mode::OPTIONAL);
    p.add("-G", "--gap_open2", "The gap open score for sequence alignment.", Mode::OPTIONAL);
    p.add("-E", "--gap_extension2", "The gap extension score for sequence alignment.", Mode::OPTIONAL);


    auto args = p.parse();

    if (!args.parsedSuccessfully()) {
        std::cerr << "Unsuccessful parse args\n";
        p.printHelpString();
        std::cout << "Exit RaMA!" << std::endl;
        return -1;
    }

    std::string data_path, output_path;
    bool save, load;
    uint_t thread_num;
    int_t match, mismatch, gap_open1, gap_open2, gap_extension1, gap_extension2;

    try {
        data_path = args["--input"];
        output_path = args["--output"];
        thread_num = args["--threads"].empty()? 0 : std::stoi(args["--threads"]);
        save = args["--save"] == "1";
        load = args["--load"] == "1";
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
    std::ios::sync_with_stdio(false);
       
	// const char* data_path = "/mnt/f/code/vs_code/RaMA/data/human.fasta";
    RareMatchPairs final_anchors;
	std::vector<SequenceInfo>* data = new std::vector<SequenceInfo>(readDataPath(data_path.c_str()));
    {
        AnchorFinder anchor_finder(*data, thread_num, output_path.c_str(), load, save);
        final_anchors = anchor_finder.lanuchAnchorSearching();
    }

    PairAligner pair_aligner(output_path, match, mismatch, gap_open1, gap_extension1, gap_open2, gap_extension2, thread_num);
    pair_aligner.alignPairSeq(*data, final_anchors);
    logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;
	delete data;
    logger.info() << "End RaMA!" << std::endl;

    return 0;

}