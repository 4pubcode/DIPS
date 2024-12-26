// #include "ds.h"
// #include "lookup_table.h"
// #include "range.h"
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "mdss.h"
#include "ppsMethod.h"
#include "util.h"
// #define BASE 4
// #define M 3
// #define RI_MAX 100
// #define II_MAX 10

// extern lookup_table *g_tables;
// extern vector<vector<CountType>> *g_bits2set;

class dPps : public ppsMethod {
   private:
    md_index *index;

   public:
    dPps(/* args */);
    ~dPps();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &results);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
};

dPps::dPps(/* args */) {}

dPps::~dPps() {}

int dPps::init(vector<CountType> elements,
                           vector<WeightType> weights) {
    index = new md_index(elements, weights);
    return 1;
}

int dPps::query(vector<CountType> &results) {
    index->sample(results);
    // vector<int> cand2;
    // for (int i=0; i<loglogn; i++) {
    //     if (dsfmt_gv_genrand_open_close() < e2p3[i]) {
    //         cand2.push_back(i);
    //     }
    // }
    // vector<int> cand1;
    // for (auto id:cand2) {
    //     WeightType q = e2p3[id];
    //     WeightType tmp = 1/(1<<id)/q;
    //     WeightType p = q>1? 1: q;
    //     int loc = geo(p);
    //     auto &arr = buckets2[id];
    //     while (loc < arr.size()) {
    //         auto e = arr[loc];
    //         if (dsfmt_gv_genrand_close_open()*q <= e2p2[e]) {
    //             cand1.push_back(e);
    //         }
    //         loc += geo(p) + 1;
    //     }
    // }
    // if (dsfmt_gv_genrand_close_open() <= rem2q) {
    //     for (auto e:rem2) {
    //         if (dsfmt_gv_genrand_close_open()*rem2q <= e2p2[e]) {
    //             cand1.push_back(e);
    //         }
    //     }
    // }

    // results.resize(0);
    // for (auto id:cand1) {
    //     WeightType q = e2p2[id];
    //     WeightType tmp = 1/(1<<id)/q;
    //     WeightType p = q>1? 1: q;
    //     int loc = geo(p);
    //     auto &arr = buckets1[id];
    //     while (loc < arr.size()) {
    //         auto e = arr[loc];
    //         if (dsfmt_gv_genrand_close_open()*q <= e2p1[e]) {
    //             results.push_back(e);
    //         }
    //         loc += geo(p) + 1;
    //     }
    // }
    // if (dsfmt_gv_genrand_close_open() <= rem1q) {
    //     for (auto e:rem1) {
    //         if (dsfmt_gv_genrand_close_open()*rem1q <= e2p1[e]) {
    //             results.push_back(e);
    //         }
    //     }
    // }
    return 1;
}

int dPps::add(WeightType addWeight, CountType idx) {
    vector<CountType> elements;
    vector<WeightType> weights;
    auto W = index->W;
    for (auto ep : index->e2p1) {
        auto e = ep.first;
        auto p = ep.second;
        elements.push_back(e);
        weights.push_back(p * W);
    }
    elements.push_back(idx);
    weights.push_back(addWeight);
    delete index;
    init(elements, weights);
    return 1;
}

int dPps::del(CountType delIdx) {
    vector<CountType> elements;
    vector<WeightType> weights;
    auto W = index->W;
    for (auto ep : index->e2p1) {
        auto e = ep.first;
        if (e != delIdx) {
            auto p = ep.second;
            elements.push_back(e);
            weights.push_back(p * W);
        }
    }
    delete index;
    init(elements, weights);
    return 1;
}

