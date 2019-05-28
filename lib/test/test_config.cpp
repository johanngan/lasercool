#include "read_config.hpp"

int main() {
    // auto params = read_config("test.cfg");
    // for(auto p: params) {
    //     std::cout << p.first << ":" << p.second << std::endl;
    // }
    double a, b, param1;
    load_params("test.cfg", {{"a", &a}, {"b", &b}, {"param1", &param1}});
    std::cout << "a:" << a << std::endl;
    std::cout << "b:" << b << std::endl;
    std::cout << "param1:" << param1 << std::endl;
}