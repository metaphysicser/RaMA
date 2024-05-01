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
#pragma once
#include "gsacak.h"
#include "kseq.h"
#include "logging.h"

#include <filesystem>
#include <random>
#include <fstream>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>


class Serializable {
public:
	virtual void serialize(std::ostream& out) const = 0;
	virtual void deserialize(std::istream& in) = 0;

	bool saveToFile(const std::string& filename) const;
	bool loadFromFile(const std::string& filename);
protected:
    template<typename T>
    void saveNumber(std::ostream& out, const T& value) const {
        out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template<typename T>
    void loadNumber(std::istream& in, T& value) {
        in.read(reinterpret_cast<char*>(&value), sizeof(T));
    }

    template<typename T>
    void saveArray(std::ostream& out, const T* array, size_t size) const {
        out.write(reinterpret_cast<const char*>(array), sizeof(T) * size);
    }

    template<typename T>
    void loadArray(std::istream& in, T* array, size_t size) {
        in.read(reinterpret_cast<char*>(array), sizeof(T) * size);
    }

    template<typename T>
    void saveVector(std::ostream& out, const std::vector<T>& vec) const {
        size_t size = vec.size();
        saveNumber(out, size); 
        for (const T& item : vec) {
            saveNumber(out, item); 
        }
    }

    template<typename T>
    void loadVector(std::istream& in, std::vector<T>& vec) {
        size_t size;
        loadNumber(in, size); 
        vec.resize(size);
        for (T& item : vec) {
            loadNumber(in, item); 
        }
    }

    template<typename T>
    void saveVector2D(std::ostream& out, const std::vector<std::vector<T>>& vec2d) const {
        size_t size = vec2d.size();
        saveNumber(out, size); 
        for (const auto& vec : vec2d) {
            saveVector(out, vec); 
        }
    }

    template<typename T>
    void loadVector2D(std::istream& in, std::vector<std::vector<T>>& vec2d) {
        size_t size;
        loadNumber(in, size); 
        vec2d.resize(size);
        for (auto& vec : vec2d) {
            loadVector(in, vec);
        }
    }
};

// Structure to hold information about a biological sequence.
struct SequenceInfo {
	std::string sequence; // The biological sequence
	std::string header;   // The header or identifier for the sequence, often from FASTA format.
	uint_t seq_len;       // The length of the sequence.

	// Constructor to initialize a SequenceInfo object with a sequence and its header.
	SequenceInfo(std::string& sequence, std::string& header);
};

// Reads sequence data from a specified file path and returns a vector of SequenceInfo objects.
// This function is typically used to read data in FASTA format.
std::vector<SequenceInfo> readDataPath(const char* ref_path, const char* query_path);

// Replaces all occurrences of the character 'N' in a given sequence string with a random nucleotide letter (A, C, G, or T).
// This is useful for dealing with unknown or ambiguous nucleotide bases in DNA sequences.
void replaceNWithRandomLetter(std::string& s);

// Ensures that a directory exists at the specified path. If the directory does not exist, it is created.
// This function is useful for setting up directories to store output files or logs.
void ensureDirExists(const std::string& path);

// Ensures that a file exists at the specified path. If the file does not exist, an empty file is created.
// This can be used to prepare output files before writing data to them.
void ensureFileExists(const std::string& path);

// Checks if a file exists at the specified path and returns true if it does, false otherwise.
// This is a utility function used to verify the existence of files before attempting to read from or write to them.
bool fileExists(const std::string& path);

// Joins two file paths, ensuring the correct path separators are used.
std::string joinPaths(const std::string& path1, const std::string& path2);

// Returns the smaller of two values.
template<typename T>
T getMinValue(const T& a, const T& b) {
    return (a < b) ? a : b; // If 'a' is less than 'b', return 'a'; otherwise, return 'b'
}

// Returns the larger of two values.
template<typename T>
T getMaxValue(const T& a, const T& b) {
    return (a > b) ? a : b; // If 'a' is greater than 'b', return 'a'; otherwise, return 'b'
}

