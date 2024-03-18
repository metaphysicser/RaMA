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

#pragma once
// modified from https://github.com/fmenozzi/argparser.git

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <iostream>
#include <cstdio>
#include <cstdlib>

// Behavior for individual args passed to add()
enum class Mode {
    REQUIRED,
    OPTIONAL,
    BOOLEAN,
};

// Object returned from parse()
class ArgMap {
private:
    std::map<std::string, std::string> args;
    bool success;

public:
    ArgMap(const std::map<std::string, std::string>& args, bool success)
        : args(args)
        , success(success) {}

    const std::string& operator[](const std::string& argstr) {
        return args[argstr];
    }

    bool parsedSuccessfully() const noexcept {
        return success;
    }
};

class Parser {
private:
    struct ArgStruct {
        std::string short_arg;
        std::string long_arg;
        std::string help_str;
        bool        bool_type;
        bool        required;
        bool        parsed;

        ArgStruct(const std::string& sa,
            const std::string& la,
            const std::string& hs,
            bool bt,
            bool rq,
            bool ps)
            : short_arg(sa)
            , long_arg(la)
            , help_str(hs)
            , bool_type(bt)
            , required(rq)
            , parsed(ps) {}
    };

    int                      m_argc;
    std::vector<std::string> m_argv;
    std::vector<ArgStruct>   m_args;

    bool m_any_adds_failed = false;

private:
    void removeEquals(std::vector<std::string>& argv) const {
        int new_argc = std::count_if(argv.begin(), argv.end(), [](const std::string& s) {
            return s.find("=") != std::string::npos;
            }) + argv.size();

            argv.reserve(new_argc);

            auto it = argv.begin();
            while (it != argv.end()) {
                auto idx = it->find("=");
                if (idx != std::string::npos) {
                    auto arg = it->substr(0, idx);
                    auto val = it->substr(idx + 1);

                    it = argv.erase(it);
                    if (!val.empty()) {
                        it = argv.insert(it, val);
                    }
                    it = argv.insert(it, arg);
                }
                if (it != argv.end()) {
                    ++it;
                }
            }
    }

    bool isMultiShortarg(const std::string& s) const noexcept {
        return s[0] == '-' && s[1] != '-' && s.size() > 2;
    }

    void expandShortargs(std::vector<std::string>& argv) const {
        int new_argc = argv.size();
        for (const auto& arg : argv) {
            if (this->isMultiShortarg(arg)) {
                new_argc += arg.size() - 2;
            }
        }

        argv.reserve(new_argc);

        auto it = argv.begin();
        while (it != argv.end()) {
            auto arg = *it;
            if (this->isMultiShortarg(arg)) {
                it = argv.erase(it);
                for (size_t i = arg.size() - 1; i > 0; i--) {
                    it = argv.insert(it, "-" + std::string(1, arg[i]));
                }
            }
            else {
                ++it;
            }
        }
    }


public:
    Parser(int argc, char* argv[]) {
        m_argv = std::vector<std::string>(argv, argv + argc);

        // Reformat argv in case --arg=val notation is used
        this->removeEquals(m_argv);

        // Expand shortargs in case -ab notation is used
        this->expandShortargs(m_argv);

        m_argc = m_argv.size();
    }

    bool add(const std::string& shortarg,
        const std::string& longarg,
        const std::string& helpstr,
        Mode m = Mode::OPTIONAL) {

        // Can't have both argstrings be empty
        if (shortarg.empty() && longarg.empty()) {
            m_any_adds_failed = true;
            return false;
        }

        // Argstrings must be formatted properly
        if (!shortarg.empty() && (shortarg.size() != 2 || shortarg[0] != '-' || shortarg[1] == '-')) {
            m_any_adds_failed = true;
            return false;
        }
        if (!longarg.empty() && (longarg.size() <= 2 || longarg[0] != '-' || longarg[1] != '-')) {
            m_any_adds_failed = true;
            return false;
        }

        // -h, --help are reserved
        if (shortarg == "-h" || longarg == "--help") {
            m_any_adds_failed = true;
            return false;
        }

        // No empty help string
        if (helpstr.empty()) {
            m_any_adds_failed = true;
            return false;
        }

        // No duplicate short/long args
        auto has_duplicate_args = [&](const ArgStruct& as) {
            return as.short_arg == shortarg || as.long_arg == longarg;
            };
        if (std::count_if(m_args.begin(), m_args.end(), has_duplicate_args) > 0) {
            m_any_adds_failed = true;
            return false;
        }

        bool booltype = (m == Mode::BOOLEAN);
        bool required = (m == Mode::REQUIRED);

        m_args.emplace_back(shortarg, longarg, helpstr, booltype, required, false);

        return true;
    }

    void printHelpString() const {
        int help_len = std::string("-h, --help").size();
        int max_len = help_len;
        int right_pad = 4;

        std::string left_pad_str = "    ";

        // Print the usage line with program name and help option
        std::cout << "Usage: " << m_argv[0] << " [-h,--help] ";
        for (const auto& as : m_args) {
            std::string sa = as.short_arg;
            std::string la = as.long_arg;

            std::string lbrak = as.required ? "" : "[";
            std::string rbrak = as.required ? "" : "]";

            if (!sa.empty()) {
                if (!la.empty()) {
                    std::cout << lbrak << sa << "," << la << rbrak << " ";
                }
                else {
                    std::cout << lbrak << sa << rbrak << " ";
                }
            }
            else {
                std::cout << lbrak << la << rbrak << " ";
            }
        }
        std::cout << "\n\n";

        // Determine the maximum argument length for formatting
        for (const auto& as : m_args) {
            int shortlen = as.short_arg.empty() ? 0 : as.short_arg.size();
            int longlen = as.long_arg.empty() ? 0 : as.long_arg.size();

            int arg_len = (shortlen && longlen) ? (shortlen + 2 + longlen) : (shortlen + longlen);
            if (arg_len > max_len)
                max_len = arg_len;
        }

        // Print the arguments section header
        std::cout << "Arguments:\n";
        std::cout << left_pad_str << "-h, --help";
        for (int i = 0; i < max_len + right_pad - help_len; i++) {
            std::cout << " ";
        }
        std::cout << "Show this help message and exit\n";

        // Print each argument with its description
        for (const auto& as : m_args) {
            std::string sa = as.short_arg;
            std::string la = as.long_arg;

            if (!sa.empty()) {
                if (!la.empty()) {
                    std::cout << left_pad_str << sa << ", " << la;
                }
                else {
                    std::cout << left_pad_str << sa;
                }
            }
            else {
                std::cout << left_pad_str << la;
            }

            int shortlen = sa.empty() ? 0 : sa.size();
            int longlen = la.empty() ? 0 : la.size();
            int arg_len = (shortlen && longlen) ? (shortlen + 2 + longlen) : (shortlen + longlen);
            for (int j = 0; j < max_len + right_pad - arg_len; j++) {
                std::cout << " ";
            }

            std::cout << as.help_str << "\n";
        }
    }


    ArgMap parse() {
        std::map<std::string, std::string> map;
        bool success = true;

        bool help_passed = false;

        if (m_any_adds_failed) {
            success = false;
        }
        else if (m_argc == 2 && (m_argv[1] == "-h" || m_argv[1] == "--help")) {
            // Check if -h, --help was passed as only arg
            this->printHelpString();
            std::exit(EXIT_SUCCESS);
        }
        else {
            // Check for rogue "="
            auto is_rogue_equal = [](const std::string& s) { return s == "="; };
            if (std::any_of(m_argv.begin(), m_argv.end(), is_rogue_equal)) {
                success = false;
            }
            else {
                // Initialize all booltype args to false and all other
                // args to the empty string
                for (const auto& arg : m_args) {
                    auto default_val = arg.bool_type ? "0" : "";

                    map[arg.short_arg] = default_val;
                    map[arg.long_arg] = default_val;
                }

                // Assign args
                for (int i = 1; i < m_argc; i++) {
                    for (auto& as : m_args) {
                        if (as.short_arg == m_argv[i] || as.long_arg == m_argv[i]) {
                            std::string val;
                            if (as.bool_type) {
                                val = "1";
                            }
                            else if (i + 1 < m_argc) {
                                val = m_argv[++i];
                            }
                            else {
                                success = false;
                            }

                            map[as.short_arg] = val;
                            map[as.long_arg] = val;

                            as.parsed = true;
                        }
                    }
                }
                map.erase("");

                if (success) {
                    // Check for required args
                    auto is_unparsed = [](const ArgStruct& as) { return as.required && !as.parsed; };
                    if (std::any_of(m_args.begin(), m_args.end(), is_unparsed)) {
                        success = false;
                    }
                }
            }
        }

        return ArgMap(map, success);
    }

    int argc() const noexcept {
        return m_argc;
    }

    const std::vector<std::string>& argv() const noexcept {
        return m_argv;
    }
};
   



