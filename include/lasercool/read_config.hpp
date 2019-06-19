#ifndef READ_CONFIG_HPP_
#define READ_CONFIG_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <unordered_map>
// Reads all name-value pairs in a config file and returns them in a map
// Default delimiter is a colon.
// Assumes comments are preceded by "#"
std::unordered_map<std::string, double> read_config(std::string,
    std::string delimiter = ":");
// Loads specified parameters into variables, passed as pointers
void load_params(std::string, std::vector< std::pair<std::string, double*> >,
    std::string delimiter = ":");

#endif