#ifndef IOTAG_HPP_
#define IOTAG_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <utility>

// Tag a file with a suffix and possibly a prefix
std::string tag_filename(std::string, std::string,
    std::string prefix = "", std::string separator = "_");
// Tag a file with multiple suffixes and possibly a prefix
std::string tag_filename(std::string, std::vector<std::string>,
    std::string prefix = "", std::string separator = "_");

// Synthesize a filename and a path
std::string fullfile(std::string, std::string);

#endif