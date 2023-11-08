/*
A vector class made for the purpose of accumulating a 
running sum per entry into the vector. 
*/

#include "aavector.hpp"

aavector::aavector(){}

void aavector::add(float val){
    aavec.push_back(val);
    sum = sum + val;
}

void aavector::print(void){
    for (std::size_t i; i < aavec.size(); i++){
        std::cout << aavec[i] << std::endl;;
    }
}

float aavector::getSum(void){
    return(sum);
}
