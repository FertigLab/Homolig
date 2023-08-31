
/*
The library function necessary for running the homolig input through cython from the python library. 
*/
#include "homolig.hpp"

const char *USAGE="usage: homolig [options] sequences.csv";
const char *HELP_INPUT="-h: help\n-i: seq_type\n-c: chains\n-m: metric\n-p: species";

std::vector<double> homolig (std::string &file, std::vector<std::string> &AA1, std::vector<std::string> &AA2){
    std::vector<double> output(AA1.size()); //this is assuming AA1.size() == AA2.size()
    aamatrix aamat(file);
    for (std::size_t i=0; i<AA1.size(); i++){
        output[i] = aamat.get_alignment(AA1[i],AA2[i]);
    }
    return(output);
}
