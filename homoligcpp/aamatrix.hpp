// aamatrix.hpp

#ifndef AAMATRIX_HPP
#define AAMATRIX_HPP

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <filesystem>
#include "aavector.hpp"


class aamatrix{
    private:
        std::vector<std::vector<float>> mat;
        std::map<std::string,int> AA_index;
    public:
        std::string AA_file;
        aamatrix(std::string file);
        // int check_mat(void);
        // void print_AA_index(void);
        float get_val(std::string AA1, std::string AA2);
        std::vector<std::string> get_kmers(std::string AA_string, int KMER_size);
        float get_alignment (std::string AA1_string, std::string AA2_string);
};

#endif
