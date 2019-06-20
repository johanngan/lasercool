#include "lasercool/readcfg.hpp"

std::unordered_map<std::string, double> read_config(std::string fname,
    std::string delimiter) {
    std::unordered_map<std::string, double> params;

    std::ifstream infile(fname);
    std::string line;
    while(std::getline(infile, line)) {
        // empty line or comment is ignored
        if(line.empty() || line.rfind("#", 0) == 0) continue;
        auto delim_pos = line.find(delimiter);
        if(delim_pos == line.npos) {
            std::cout 
                << "Error in read_config(): encountered line without delimiter."
                << std::endl;
            return params;
        }
        std::string key = line.substr(0, delim_pos);
        double value = std::stod(line.substr(delim_pos+1));
        params[key] = value;
    }
    infile.close();
    return params;
}

void load_params(std::string fname,
    std::vector< std::pair<std::string, double*> > params,
    std::string delimiter) {
    auto all_params = read_config(fname, delimiter);
    for(auto p: params) {
        *p.second = all_params[p.first];
    }
}