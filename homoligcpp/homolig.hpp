#ifndef HOMOLIG_HPP
#define HOMOLIG_HPP

#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include "aamatrix.hpp"

extern const char *USAGE;
extern const char *HELP_INPUT;

std::vector <double> homolig(std::string &file,std::vector<std::string> &AA1, std::vector<std::string> &AA2);
// void homolig(std::string file);

#endif
