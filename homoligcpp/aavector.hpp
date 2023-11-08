#ifndef AAVECTOR_HPP
#define AAVECTOR_HPP

#include <vector>
#include <iostream>

class aavector{
    private:
        std::vector<float> aavec;
        float sum=0;
    public:
        aavector(void);
        void add(float val);
        float getSum(void);
        void print(void);
};

#endif
