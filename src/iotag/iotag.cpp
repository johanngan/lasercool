#include "lasercool/iotag.hpp"

std::string tag_filename(std::string filename, std::string suffix,
    std::string prefix, std::string separator) {
    std::string tagged = filename.insert(filename.rfind("."),
        separator + suffix);
    if(!prefix.empty()) {
        tagged = prefix + separator + tagged;
    }
    return tagged;
}

std::string tag_filename(std::string filename,
    std::vector<std::string> suffixes, std::string prefix,
    std::string separator) {
    std::string tagged = filename;
    for(auto s: suffixes) {
        tagged = tag_filename(tagged, s, "", separator);
    }
    if(!prefix.empty()) {
        tagged = prefix + separator + tagged;
    }
    return tagged;
}

std::string fullfile(std::string filename, std::string dir) {
    return dir + "/" + filename;
}

std::pair<std::string, std::string> fileparts(std::string path) {
    std::string dir = path.substr(0, path.rfind("/"));
    std::string filename = path.substr(path.rfind("/") + 1);
    return std::make_pair(filename, dir);
}