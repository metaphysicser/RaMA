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
#include "utils.h"
// If the C4996 error occurs in Visual Studio, disabling this error in the IDE will allow successful execution.
KSEQ_INIT(int, read)

bool Serializable::saveToFile(const std::string& filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        logger.error() << "Cannot open file for writing: " + filename << std::endl;
        return false;
    }
    serialize(out);
    out.close();
    return true;
}

bool Serializable::loadFromFile(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        logger.error() << "Cannot open file for reading: " + filename << std::endl;
        return false;
    }
    deserialize(in);
    in.close();
    return true;
}

// SequenceInfo constructor: Initializes a sequence information instance with a sequence and its header.
SequenceInfo::SequenceInfo(std::string& seq_value, std::string& header_value) {
    sequence = seq_value; // Assign the sequence value.
    header = header_value; // Assign the header value.
    seq_len = sequence.length(); // Calculate and store the length of the sequence.
}

// Function to read a single sequence from a file
SequenceInfo readSingleSequence(const char* filepath) {
    FILE* f_pointer = fopen(filepath, "r");
    if (!f_pointer) {
        logger.error() << "Failed to open the fasta file at: " << filepath << std::endl;
        logger.error() << "RaMA Exit!" << std::endl;
        exit(EXIT_FAILURE);
    }

    kseq_t* seq = kseq_init(fileno(f_pointer));
    if (kseq_read(seq) < 0) {
        logger.error() << "Failed to read a sequence from: " << filepath << std::endl;
        logger.error() << "RaMA Exit!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string seq_data = seq->seq.s;
    replaceNWithRandomLetter(seq_data);
    std::string seq_name = seq->name.s;
    if (seq->comment.l) seq_name += " " + std::string(seq->comment.s);

    kseq_destroy(seq);
    fclose(f_pointer);

    return SequenceInfo(seq_data, seq_name);
}

// Reads sequence data from a specified file path and constructs a vector of SequenceInfo objects.
std::vector<SequenceInfo> readDataPath(const char* ref_path, const char* query_path) {
    // Check if the specified file exists.
    if (!fileExists(ref_path)) {
        logger.error() << "Cannot find the reference fasta data, please check the data path." << std::endl;
        logger.error() << "RaMA Exit!" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!fileExists(query_path)) {
        logger.error() << "Cannot find the query fasta data, please check the data path." << std::endl;
        logger.error() << "RaMA Exit!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<SequenceInfo> data;
    data.push_back(readSingleSequence(ref_path));
    data.push_back(readSingleSequence(query_path));

    // Ensure there are exactly two sequences in the file.
    if (data.size() != 2) {
        logger.error() << "The fasta file should only have 2 sequences, but there are "
            << data.size() << " in it." << std::endl;
        logger.error() << "RaMA Exit!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Log the lengths of the first and second sequences.
    logger.info() << "The first sequence length is " << data[0].seq_len << std::endl;
    logger.info() << "The second sequence length is " << data[1].seq_len << std::endl;

    return data; // Return the vector of SequenceInfo objects.
}

// Replaces 'n' characters in a given string with random nucleotide letters.
void replaceNWithRandomLetter(std::string& s) {
    if (s.empty()) return; // Do nothing if the string is empty.

    std::random_device rd; // Obtain a random number from hardware.
    std::mt19937 gen(rd()); // Seed the generator.
    std::uniform_int_distribution<> distr(0, 3); // Define the range for nucleotide indices.
    const char* bases = "ACGT"; // Nucleotide bases.

    // Iterate through each character of the string.
    for (auto& ch : s) {
        if (ch == 'n' || ch == 'N') {
            ch = bases[distr(gen)]; // Replace 'n' with a random nucleotide.
        }
        ch = std::toupper(static_cast<unsigned char>(ch));

    }
}

void ensureDirExists(const std::string& path) {
    std::filesystem::path fsPath(path); 

    if (!std::filesystem::exists(fsPath)) {
        try {
            if (std::filesystem::create_directories(fsPath)) {
            }
            else {
                std::cerr << "Failed to create directory: " << fsPath << std::endl;
            }
}
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
        }
    }
    else if (!std::filesystem::is_directory(fsPath)) {
        std::cerr << "Path exists but is not a directory: " << fsPath << std::endl;
    }
}

// Ensures that a file exists at the specified path, creating an empty file if necessary.
void ensureFileExists(const std::string& path) {
    std::filesystem::path fsPath(path);
    // Check if the file does not exist.
    if (!std::filesystem::exists(fsPath)) {
        // Create a new empty file.
        std::ofstream create_file(fsPath);
        // Ensure the file was successfully created.
        assert(create_file.is_open());
    }
}

// Checks if a file exists at the specified path.
bool fileExists(const std::string& path) {
    std::filesystem::path fsPath(path);
    // Check if the path exists and is not a directory.
    return std::filesystem::exists(fsPath) && !std::filesystem::is_directory(fsPath);
}

// Joins two paths into a single path, handling cases with or without trailing/leading separators.
std::string joinPaths(const std::string& path1, const std::string& path2) {
    // Use std::filesystem::path to handle the concatenation and normalization of paths.
    std::filesystem::path result = std::filesystem::path(path1) / std::filesystem::path(path2);
    // Convert the resulting path back to a string.
    return result.string();
}
