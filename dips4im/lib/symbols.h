#pragma once
#include <iostream>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>

#include "Argument.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "SFMT/dSFMT/dSFMT.h"
#include "alias.h"
#define MY_C (1)

using namespace std;

typedef unsigned int uint;
typedef uint CountType;
typedef double WeightType;
typedef int NameType;
typedef pair<CountType, WeightType> element;
const CountType PRESERVE_SUBSET = 3;
const CountType INTBITS = 32;
const double numSubsetSizeBoundMul = 1.5;
uint CountTypeMAX = UINT32_MAX;
const double weightSubsetGap = 4;

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

// #include "bss_graph.h"
// #include "dss_graph.h"
#include "graph.h"
// #include "hss_graph.h"
#include "logger.h"
// #include "nss_graph.h"
#include "timer.h"
#include "tools.h"
// #include "bringmannSubsetSampling.h"
// #include "hybridSubsetSampling.h"
// #include "naiveSubsetSampling.h"
// #include "dynamicSubsetSampling.h"
#include "alg.h"
#include "commonFunc.h"
#include "commonStruct.h"
#include "hyperGraph.h"