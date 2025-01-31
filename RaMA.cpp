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

	p.add("-r", "--reference", "Reference FASTA file path containing the reference sequences for alignment.", Mode::REQUIRED);
	p.add("-q", "--query", "Query FASTA file path containing the query sequences for alignment.", Mode::REQUIRED);
	p.add("-o", "--output", "Output directory path for saving alignment results and additional files.", Mode::REQUIRED);

	p.add("-t", "--threads", "Number of threads for the alignment process. Defaults to the number of available cores if unspecified.", Mode::OPTIONAL);

	p.add("-s", "--save", "Saves anchor binary files to the output directory for future use, including SA, LCP, and Linear Sparse Table.", Mode::BOOLEAN);
	p.add("-l", "--load", "Loads existing anchor binary files from the output directory to skip SA, LCP, and Linear Sparse Table construction.", Mode::BOOLEAN);

	p.add("-c", "--max_match_count", "Maximum number of rare matches to use for anchor finding. Altering this value is generally not recommended.", Mode::OPTIONAL);

	p.add("-m", "--match", "Match score for sequence alignment. Lower values favor matching characters. Default is 0.", Mode::OPTIONAL);
	p.add("-x", "--mismatch", "Mismatch penalty. Higher values penalize mismatches more. Default is 3.", Mode::OPTIONAL);
	p.add("-g", "--gap_open1", "Penalty for initiating a short gap. Key for handling different gap lengths. Default is 4.", Mode::OPTIONAL);
	p.add("-e", "--gap_extension1", "Penalty for extending a short gap. Less severe than gap opening penalty. Default is 2.", Mode::OPTIONAL);
	p.add("-G", "--gap_open2", "Penalty for initiating a long gap. Aims to manage long gaps strategically. Default is 12.", Mode::OPTIONAL);
	p.add("-E", "--gap_extension2", "Penalty for extending a long gap. Provides a lenient approach to long gap management. Default is 1.", Mode::OPTIONAL);

	p.add("-a", "--sam_output", "Whether to output in SAM format. If specified, results will be saved in SAM format.", Mode::BOOLEAN);
	p.add("-p", "--paf_output", "Whether to output in PAF format. If specified, results will be saved in PAF format.", Mode::BOOLEAN);

	// Parse the command line arguments using the Parser instance
	auto args = p.parse();

	// Check if the parsing was successful
	if (!args.parsedSuccessfully()) {
		std::cerr << "Unsuccessful parse args\n";
		p.printHelpString(); // Print the help string if parsing was unsuccessful
		std::cout << "Exit RaMA!" << std::endl;
		return -1; // Exit with an error code
	}

	// Initialize variables for storing command line arguments
	std::string ref_path, query_path, output_path;
	bool save, load, sam_output, paf_output;
	uint_t thread_num, max_match_count;
	int_t match, mismatch, gap_open1, gap_open2, gap_extension1, gap_extension2;

	try {
		// Assign the parsed values to variables, with defaults where necessary
		ref_path = args["--reference"];
		query_path = args["--query"];
		output_path = args["--output"];
		thread_num = args["--threads"].empty() ? std::thread::hardware_concurrency() : std::stoi(args["--threads"]);
		save = args["--save"] == "1";
		load = args["--load"] == "1";
		sam_output = args["--sam_output"] == "1";
		paf_output = args["--paf_output"] == "1";
		max_match_count = getMaxValue(args["--max_match_count"].empty() ? 100 : std::stoi(args["--max_match_count"]), 2);
		match = args["--match"].empty() ? 0 : std::stoi(args["--match"]);
		mismatch = args["--mismatch"].empty() ? 3 : std::stoi(args["--mismatch"]);
		gap_open1 = args["--gap_open1"].empty() ? 4 : std::stoi(args["--gap_open1"]);
		gap_extension1 = args["--gap_extension1"].empty() ? 2 : std::stoi(args["--gap_extension1"]);
		gap_open2 = args["--gap_open2"].empty() ? 12 : std::stoi(args["--gap_open2"]);
		gap_extension2 = args["--gap_extension2"].empty() ? 1 : std::stoi(args["--gap_extension2"]);
	}
	catch (std::exception& e) {
		// Catch and report any errors during argument processing
		std::cerr << "Error: " << e.what() << std::endl;
		std::cerr << "Unsuccessful parse args\n";
		p.printHelpString();
		std::cout << "Exit RaMA!" << std::endl;
		return -1;
	}

	// Set the output directory for logging
	logger.setDir(output_path);
	logger.info() << "Start RaMA!" << std::endl;

	RareMatchPairs final_anchors;
	// Load sequences from the input data path
	std::vector<SequenceInfo>* data = new std::vector<SequenceInfo>(readDataPath(ref_path.c_str(), query_path.c_str()));
	{
		// Initialize AnchorFinder with the provided arguments and find anchors
		AnchorFinder anchor_finder(*data, output_path.c_str(), thread_num, load, save, max_match_count);
		final_anchors = anchor_finder.lanuchAnchorSearching();
	}
	// final_anchors.clear();
	// std::cout << final_anchors.size() << std::endl;
	// Initialize PairAligner with the parsed arguments and align the sequences
	PairAligner pair_aligner(output_path, match, mismatch, gap_open1, gap_extension1, gap_open2, gap_extension2, thread_num);
	pair_aligner.alignPairSeq(*data, final_anchors, sam_output, paf_output);

	// Log the maximum memory used during the process
	logger.info() << "Max memory used is " << logger.getMaxMemoryUsed() << std::endl;

	// Clean up allocated memory for data
	delete data;

	// Log the completion of the process
	logger.info() << "End RaMA!" << std::endl;

	// Return success
	return 0;
}