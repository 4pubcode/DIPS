#pragma once
// #include "symbols.h"
class ppsMethod {
   private:
    /* data */
   public:
    ppsMethod(/* args */);
    ~ppsMethod();
    virtual int init(vector<CountType> elements, vector<WeightType> weights) {return 0;}
    virtual int query(vector<CountType>& results) {return 0;}
    virtual int add(WeightType& addWeight, CountType idx) {return 0;}
    virtual int del(CountType& delIdx) {return 0;}
};

ppsMethod::ppsMethod(/* args */) {}

ppsMethod::~ppsMethod() {}
