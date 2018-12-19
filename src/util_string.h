#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <regex>
#include <set>

using namespace std;

void reverseSequence(std::string &seq);

void reverseComplement(std::string &seq);

std::vector<std::string> split_string(const std::string &s,
                                      const std::string &delim);
pair<int, string> strRepeatEasy(const string &str);

int countCharInStr(string str, char i);
string condense_repeat_string(string expanded_string);

int get_barcode_index_bc(string bc, const std::vector<std::string> &barcodes_list);

int get_barcode_index_bc_exact(string bc, const std::vector<std::string> &barcodes_list);

void read_barcode_list(string fn, std::vector<std::string> &barcodes_list);

string trim_barcode(string barcode);