#pragma once
#include <iostream>
using namespace std;
#define MY_C (1)
#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <string>

#include "SFMT/dSFMT/dSFMT.c"
#include "SFMT/dSFMT/dSFMT.h"
// typedef unsigned int uint;
// typedef uint CountType;
typedef uint32_t CountType;
typedef double WeightType;
typedef uint32_t NameType;
typedef pair<NameType, WeightType> element;
const CountType PRESERVE_SUBSET = 3;
const CountType INTBITS = 32;
const double numSubsetSizeBoundMul = 1.5;
const double weightSubsetGap = 4;
CountType CountTypeMAX = UINT32_MAX;
#include <set>
#include <unordered_map>

#include "logger.h"
#include "timer.h"
#include "tools.h"
typedef struct {
    element *addr;
    WeightType maxw;
    CountType size;
    CountType lastIdx;
} subsetInfo;
struct node {
    CountType id;
    double w;
    node() {}
    node(CountType _id, double _w) : id(_id), w(_w) {}
    node(const node &tmp) {
        id = tmp.id;
        w = tmp.w;
    }
    node &operator=(const node &tmp) {
        id = tmp.id;
        w = tmp.w;
        return *this;
    }
};

#include "Argument.h"
