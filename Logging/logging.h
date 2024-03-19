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
#pragma once

#include "gsacak.h"
#include "utils.h"

#include <cassert>
#include <stdint.h>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <fstream>
#include <sys/resource.h>
#include <unistd.h>

// PerformanceMonitor class is used for monitoring and retrieving performance metrics,
// such as memory usage and execution time.
class PerformanceMonitor {
public:
    // Constructor initializes the monitoring process.
    PerformanceMonitor();

    // Returns a string containing various performance metrics.
    std::string getPerformanceMetrics();

    // Retrieves the maximum memory usage recorded during the lifetime of the monitor.
    std::string getMaxMemory() const;

    // Resets the monitoring data to start fresh measurements.
    void reset();

private:
    // Formats the memory usage from bytes to a more readable string format.
    std::string formatMemoryUsage(uint_t bytes) const;

    // Stores the start time of the performance monitoring for calculating elapsed time.
    std::chrono::steady_clock::time_point start_time;

    // Records the maximum memory usage observed.
    uint_t max_memory;
};

// LogLevel enum defines different levels of logging, such as error, info, and debug.
enum LogLevel { error, info, debug };

// Logger class is designed for logging messages with various levels of severity.
// It extends std::streambuf and std::ostream for custom logging behavior.
class Logger : public std::streambuf, public std::ostream {
private:

    std::string program_name;

    // Directory where log files are stored.
    std::string dir;

    // Path to the current log file.
    std::string log_file;

    // Directory where backups of log files are stored.
    std::string backup_dir;

    // Indicates whether to also print log messages to the standard output.
    bool add_cout;

    // Output stream for writing log messages to a file.
    std::ofstream* os;

    // PerformanceMonitor instance for recording and logging performance metrics.
    PerformanceMonitor* monitor;

    // Current logging level of a message being processed.
    LogLevel cur_level;

    // Maximum logging level to output; messages above this level are ignored.
    LogLevel max_level;

public:
    // Constructor initializes the Logger with a directory for logs, program name,
    // whether to add console output, and the maximum log level.
    explicit Logger(const std::string& _program_name, bool _add_cout = true, LogLevel _max_level = LogLevel::info);

    // Copy constructor is deleted to prevent copying of Logger instances.
    Logger(const Logger&) = delete;

    void setDir(std::string& dir);

    // Adds a new log entry, potentially rotating the log file if necessary.
    void addNewLog();

    // Backs up the current log file to a specified backup directory.
    void backup();

    // Override from std::streambuf, handles the buffering of log messages.
    int overflow(int c) override;

    // Forces flushing the buffered log messages to the output stream.
    void forceFlush();

    // Sets the current log level to error, info, or debug, and returns a reference to the logger.
    Logger& error();
    Logger& info();
    Logger& debug();

    bool isDebugEnabled();

    std::string getMaxMemoryUsed() const;

    // Destructor cleans up resources, such as closing file streams.
    ~Logger() override;
};

// Global logger instance for use throughout the application.
extern Logger logger;
