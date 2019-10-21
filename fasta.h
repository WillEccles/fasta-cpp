/*
 * fasta.h
 *
 * Author: Will Eccles
 * Date: 2019-10-17
 * Description: Implements a simple header-only FASTA file parser.
 * 
 * See also:
 *   https://en.wikipedia.org/wiki/FASTA_format
 *
 * Copyright 2019 Will Eccles
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

#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <iostream>

// Macro to convert a character to uppercase
#define TO_UPPER(c) (((c) >= 'a' && (c) <= 'z') ? ((c) - ('a' - 'A')) : (c))

// Valid codes. Source: Wikipedia
const char* VALID_CODES = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*-";

class FASTAFile {
    public:
        FASTAFile(): file("") {}
        FASTAFile(const std::string& filename): file(filename) {
            if (!open(filename)) {
                throw std::runtime_error("Error opening file: " + filename + "!");
            }
        }

        ~FASTAFile() { infile.close(); }

        // returns false if opening the file failed
        // use this only if you used the default constructor.
        bool open(const std::string& filename) {
            file = filename;
            infile = std::ifstream(file);
            if (!infile) return false;
            std::string tmp;
            std::string head = "";
            int lc = 0;
            // deal with comments at the top of the file
            while (std::getline(infile, tmp)) {
                if (tmp[0] != '>' && tmp[1] != ';') {
                    // line length
                    line_nt = tmp.size();
                    break;
                } else {
                    head += tmp;
                    lc++;
                }
            }
            header_len = head.size() + lc; // +1 for newline character
            return true;
        }

        // closes the file
        void close() {
            infile.close();
        }

        // gets a string of nucleotides from start to end, inclusive;
        // so specifying 1, 2 would get 2nt
        // 'caps' defaults to false. If specified, this uppercases all nucleotides
        std::string get_sequence(std::size_t start, std::size_t end, bool caps = false) {
            std::string ret;
            char tmp;
            infile.seekg(seq_start(start));
            std::size_t count = 0;
            while (count < (end - start) + 1) {
                tmp = infile.get();

                if (tmp == std::ifstream::traits_type::eof()) {
                    throw std::runtime_error("End coordinate out of bounds");
                }

                if (nullptr != std::strchr(VALID_CODES, tmp)) {
                    if (caps) {
                        ret += TO_UPPER(tmp);
                    } else {
                        ret += tmp;
                    }
                    count++;
                }
            }
            return ret;
        }

    private:
        std::size_t header_len;
        std::size_t line_nt;
        std::string file;
        std::ifstream infile;

        // get the *actual* starting byte to seek to
        inline std::size_t seq_start(std::size_t start) {
            return (header_len + ((start / line_nt) * (line_nt + 1)) + (start % line_nt)) - 1;
        }
};

#endif /* FASTA_H */
