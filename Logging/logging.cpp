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
 // Created: 2024-01-30
#include "logging.h"

// Constructor: Initializes the performance monitor by resetting its metrics.
PerformanceMonitor::PerformanceMonitor() {
    reset();
}

// Retrieves performance metrics, including elapsed time and memory usage.
std::string PerformanceMonitor::getPerformanceMetrics() {
    // Capture the current time to calculate elapsed time.
    auto endTime = std::chrono::steady_clock::now();
    // Calculate elapsed time in seconds since the monitor was reset or initialized.
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(endTime - start_time).count();

    std::stringstream ss;
    // Format elapsed time into HH:MM:SS format for readability.
    ss << std::setw(2) << std::setfill('0') << elapsed / 3600 << ":"
        << std::setw(2) << (elapsed / 60) % 60 << ":"
        << std::setw(2) << elapsed % 60 << " ";

    // Platform-specific memory usage retrieval.
#if defined(_WIN32) || defined(_WIN64)
    PROCESS_MEMORY_COUNTERS memCounter;
    // Attempt to get the current process's memory usage on Windows.
    if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof(memCounter))) {
        // Format and append memory usage to the string stream.
        ss << formatMemoryUsage(memCounter.WorkingSetSize);
        // Update max_memory if the current usage is greater.
        if (memCounter.WorkingSetSize > max_memory) {
            max_memory = memCounter.WorkingSetSize;
        }
    }
    else {
        ss << "Memory usage unavailable"; // Handle case where memory info cannot be retrieved.
    }
#else
    struct rusage usage;
    // Attempt to get the current process's memory usage on Unix/Linux.
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // Format and append memory usage (convert kilobytes to bytes first).
        ss << formatMemoryUsage(usage.ru_maxrss * 1024L);
        // Update max_memory if the current usage is greater.
        if (usage.ru_maxrss * 1024L > max_memory) {
            max_memory = usage.ru_maxrss * 1024L;
        }
    }
    else {
        ss << "Memory usage unavailable"; // Handle case where memory info cannot be retrieved.
    }
#endif

    return ss.str(); // Return the formatted string.
    }

// Formats memory usage from bytes to a human-readable string with appropriate units.
std::string PerformanceMonitor::formatMemoryUsage(uint_t bytes) const {
    std::stringstream ss;
    // Choose appropriate units (B, KB, MB, GB) based on the size and format the string.
    if (bytes < 1024) {
        ss << bytes << " B";
    }
    else if (bytes < 1048576) {
        ss << std::fixed << std::setprecision(1) << bytes / 1024.0 << " KB";
    }
    else if (bytes < 1073741824) {
        ss << std::fixed << std::setprecision(1) << bytes / 1048576.0 << " MB";
    }
    else {
        ss << std::fixed << std::setprecision(1) << bytes / 1073741824.0 << " GB";
    }
    return ss.str();
}

// Resets the performance monitor to start a new monitoring session.
void PerformanceMonitor::reset() {
    // Reset start time to the current time.
    start_time = std::chrono::steady_clock::now();
    // Reset maximum memory usage observed to 0.
    max_memory = 0;
}

// Retrieves the maximum memory usage observed in a formatted string.
std::string PerformanceMonitor::getMaxMemory() const {
    std::stringstream ss;
    ss << "MaxMemory: ";
    // Format the maximum memory usage into a readable string with appropriate units.
    if (max_memory < 1024) {
        ss << max_memory << " B";
    }
    else if (max_memory < 1048576) {
        ss << std::fixed << std::setprecision(1) << max_memory / 1024.0 << " KB";
    }
    else if (max_memory < 1073741824) {
        ss << std::fixed << std::setprecision(1) << max_memory / 1048576.0 << " MB";
    }
    else {
        ss << std::fixed << std::setprecision(1) << max_memory / 1073741824.0 << " GB";
    }
    return ss.str(); // Return the formatted maximum memory usage.
}


// Logger constructor: Initializes the logger with specified directory, program name, 
// whether to add console output, and the maximum log level.
Logger::Logger(const std::string& _program_name, bool _add_cout, LogLevel _max_level) :
    std::ostream(this), // Initialize the ostream part of this logger to use this instance as the buffer.
    program_name(_program_name), // Set the program name for constructing log file names.
    add_cout(_add_cout), // Set whether to also output log messages to the console.
    os(new std::ofstream()), // Allocate a new ofstream for log file writing.
    monitor(new PerformanceMonitor()), // Create a new PerformanceMonitor instance for tracking performance metrics.
    cur_level(LogLevel::info), // Set the current log level to INFO by default.
    max_level(_max_level) { // Set the maximum log level for filtering messages.
}

void Logger::setDir(std::string& dir) {
    this->dir = dir;
    this->log_file = joinPaths(dir, (program_name + ".log"));  // Construct the log file name using the program name.
    this->backup_dir = joinPaths(dir, "old_logs/"); // Set the directory for storing backup log files.
    ensureDirExists(dir); // Ensure the log directory exists.
    ensureDirExists(backup_dir); // Ensure the backup directory exists.
    addNewLog(); // Add a new log entry (and potentially create a new log file).
}

// addNewLog: Closes the current log file if open and backs it up if it exists before creating a new log file.
void Logger::addNewLog() {
    if (os && os->is_open()) // If there's an open log file,
        os->close(); // close it.
    if (fileExists(log_file)) // If the log file already exists,
        backup(); // back it up.

    os->open(log_file, std::ios::out | std::ios::app); // Open a new or existing log file for appending.
    assert(os && os->is_open()); // Assert that the file is successfully opened.
}

// backup: Creates a backup of the current log file with a timestamp and potentially a random number to avoid overwrites.
void Logger::backup() {
    if (os && os->is_open()) { // If there's an open log file,
        os->close(); // close it.
    }
    ensureDirExists(backup_dir); // Ensure the backup directory exists.

    // Generate a backup file name using the current date and time.
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t_now), "%Y%m%d_%H%M%S");
    std::string timestamp = ss.str();
    std::string newBackupFile = backup_dir + "/" + timestamp + "_" + std::to_string(rand()) + ".log";

    // Perform the backup by copying the log file to the backup file.
    std::ifstream src(log_file, std::ios::binary);
    std::ofstream dst(newBackupFile, std::ios::binary);
    if (!src.is_open() || !dst.is_open()) {
        std::cerr << "Failed to open source or destination file for backup.\n";
        return;
    }
    dst << src.rdbuf();

    src.close(); // Close the source file.
    dst.close(); // Close the destination file.

    // Remove the original log file after backup.
    if (std::remove(log_file.c_str()) != 0) {
        std::cerr << "Failed to remove original log file: " << log_file << std::endl;
    }
}

// overflow: Custom behavior for handling character overflow. This function is called by ostream 
// when it needs to output a character. It routes the character to the appropriate destination(s).
int Logger::overflow(int c) {
    if (add_cout && cur_level == LogLevel::info)
        std::cout << char(c); // Output to console if enabled and current level is info.
    if (cur_level == LogLevel::error)
        std::cerr << char(c); // Output to standard error if current level is error.
    if (os->is_open() && cur_level <= max_level)
        *os << char(c); // Write to the log file if it's open and the current level is within the max level.
    if (c == '\n') {
        forceFlush(); // Flush the streams after a newline character.
    }
    return 0;
}

// forceFlush: Forces flushing of the output streams to ensure all messages are written.
void Logger::forceFlush() {
    if (cur_level <= LogLevel::info)
        std::cout.flush(); // Flush console output.
    os->flush(); // Flush file output.
}

// info/debug/error: These functions set the current logging level and prepend the log message with appropriate level indicator and performance metrics.
Logger& Logger::info() {
    cur_level = LogLevel::info;
    *this << monitor->getPerformanceMetrics() << " INFO: ";
    return *this;
}

Logger& Logger::error() {
    cur_level = LogLevel::error;
    *this << monitor->getPerformanceMetrics() << " ERROR: ";
    return *this;
}

Logger& Logger::debug() {
    cur_level = LogLevel::debug;
    *this << monitor->getPerformanceMetrics() << " DEBUG: ";
    return *this;
}

bool Logger::isDebugEnabled() {
    return this->max_level >= LogLevel::debug;
}

std::string Logger::getMaxMemoryUsed() const{
    return monitor->getMaxMemory();
}

// Destructor: Cleans up resources by closing the log file and deleting dynamically allocated objects.
Logger::~Logger() {
    if (os) {
        if (os->is_open()) {
            os->close(); // Close the log file if open.
        }
        delete os; // Delete the ofstream object.
        os = nullptr;
    }
    if (monitor) {
        delete monitor; // Delete the PerformanceMonitor object.
    }
}
