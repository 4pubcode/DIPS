// #include "ds.h"
// #include "lookup_table.h"
// #include "range.h"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "ppsMethod.h"
#include "util.h"
#include "mdss.h"
// #define BASE 4
// #define M 3
// #define RI_MAX 100
// #define II_MAX 10

// extern lookup_table *g_tables;
// extern vector<vector<CountType>> *g_bits2set;

class dPps : public ppsMethod {
   private:
    // vector<vector<CountType>> buckets1;
    // WeightType rem1q;
    // vector<CountType> rem1;
    // unordered_map<CountType, WeightType> e2p1;
    // vector<vector<CountType>> buckets2;
    // vector<CountType> rem2;
    // WeightType rem2q;
    // vector<WeightType> e2p2;
    // vector<CountType> e3;
    // vector<WeightType> e2p3;
    // lookup_table *tables;
    // vector<vector<int>> *bits2set;
    // WeightType W;
    md_index *index;
    // int n, logn, loglogn;
    

   public:
    dPps(/* args */);
    ~dPps();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &results);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
};

dPps::dPps(/* args */){}

dPps::~dPps() {
}

int dPps::init(vector<CountType> elements,
                           vector<WeightType> weights) {
    index = new md_index(elements, weights);
    // n = elements.size();
    // if (n == 0) return 0;
    // logn = static_cast<int>(ceil(log2(n)));
    // loglogn = logn>0 ? static_cast<int>(ceil(log2(logn))) : 0;
    // buckets1.resize(logn);
    // buckets2.resize(loglogn);
    // e2p2.resize(logn);
    // e3.resize(loglogn);
    // e2p3.resize(loglogn);

    // W= 0;
    // for (int i=0; i<n; i++) {
    //     W+=weights[i];
    // }

    // for (int i=0; i<logn; i++) {
    //     e2p2[i] = 0;
    // }
    // rem1q = 0;
    // for (int i=0 ;i<n;i++) {
    //     auto e = elements[i];
    //     auto p = weights[i]/W;
    //     e2p1[e] = p;
    //     size_t id = -ceil(log2(p));
    //     if (id<logn) {
    //         buckets1[id].push_back(e);
    //         e2p2[id] = 1-(1-e2p2[id])*(1-p);
    //     } else {
    //         rem1.push_back(e);
    //         rem1q = 1-(1-rem1q)*(1-p);
    //     }
    // }

    // for (int i=0; i<loglogn; i++) {
    //     e2p2[i] = 0;
    // }
    // rem2q = 0;
    // for (int i=0; i<logn; i++) {
    //     auto e = i;
    //     auto p = e2p2[i];
    //     size_t id = -ceil(log2(p));
    //     if (id<loglogn) {
    //         buckets2[id].push_back(e);
    //         e2p3[id] = 1-(1-e2p3[id])*(1-p);
    //     } else {
    //         rem2.push_back(e);
    //         rem2q = 1-(1-rem2q)*(1-p);
    //     }
    // }
    // for (int i=0; i<loglogn; i++) {
    //     if (e2p3[i]!=0) {
    //         e3.push_back(i);
    //     }
    // }
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
    for (auto ep:index->e2p1) {
        auto e = ep.first;
        auto p = ep.second;
        elements.push_back(e);
        weights.push_back(p*W);
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
    for (auto ep:index->e2p1) {
        auto e = ep.first;
        if (e!=delIdx) {
            auto p = ep.second;
            elements.push_back(e);
            weights.push_back(p*W);
        }
    }
    delete index;
    init(elements, weights);
    return 1;
}


// #include "lookup_table.h"
// #include "subsetSamplingMethod.h"
// #define BASE 4
// #define M 3
// #define RI_MAX 100
// #define II_MAX 10

// extern lookup_table *g_tables;
// extern vector<vector<CountType>> *g_bits2set;

// class mdSubsetSampling : public subsetSamplingMethod {
//    private:
//     CountType n;
//     CountType size;
//     pair<CountType, WeightType> *nameProb;
//     CountType *dict;
//     vector<CountType> *ranges;
//     WeightType *wr;
//     mdSubsetSampling **intervals;
//     WeightType *md;

//     // lookup_table *tables;
//     // vector<vector<int>> *bits2set;
//     WeightType W;
//     int f;

//    public:
//     mdSubsetSampling(/* args */);
//     ~mdSubsetSampling();
//     int init(vector<CountType> elements, vector<WeightType> weights);
//     int query(vector<CountType> &results);
//     int add(WeightType addWeight, CountType idx);
//     int del(CountType delIdx);
// };

// int mdSubsetSampling::init(vector<CountType> elements,
//                            vector<WeightType> weights) {
//     auto n = elements.size();
//     if (n == 0) return 0;
//     if (n > CountTypeMAX / 2)
//         size = CountTypeMAX;
//     else
//         size = n * 2;
//     nameProb = new pair<CountType, WeightType>[size];
//     dict = new CountType[size];
//     ranges = new vector<CountType>[RI_MAX];
//     wr = new WeightType[RI_MAX];
//     for (CountType i=0; i<RI_MAX; i++) {
//         wr[i] = 0;
//     }
//     intervals = new mdSubsetSampling*[II_MAX];
//     md = new WeightType[II_MAX];
//     for (CountType i=0; i<II_MAX; i++) {
//         intervals[i] = nullptr;
//         md[i] = 0;
//     }
//     for (CountType i = 0; i < n; i++) {
//         nameProb[i] = make_pair(elements[i], weights[i]);
//         dict[elements[i]] = i;
//     }
//     if (n < M) {
//         if (g_tables == nullptr) {
//             g_tables = new lookup_table(M - 1, BASE);
//             for (int i = 0; i < 256; i++) {
//                 auto tmp = i;
//                 for (CountType j = 0; j < 8; j++) {
//                     if (tmp % 2 == 1) {
//                         (*g_bits2set)[i].push_back(j);
//                     }
//                     tmp >>= 1;
//                 }
//             }
//         }
//         // tables = g_tables;
//         // bits2set = g_bits2set;
//         f = encode(weights, g_tables->A);
//         return 1;
//     }
//     // tables = nullptr;
//     // bits2set = nullptr;
//     f = -1;
//     int logn = static_cast<int>(ceil(mylog(n)));
//     W = 0;
//     for (CountType i=0; i<n; i++) {
//         W+=weights[i];
//         auto Ri = getRi(weights[i]);
//         auto Ii = getIi(Ri, logn);
//         if (ranges[Ri].empty())
//         if (intervals[Ii] == nullptr) {
//             intervals[Ii]
//         }
//     }
// }
// int mdSubsetSampling::query(vector<CountType> &results) {
//     if (n < M) {
//         auto idxs = (*g_bits2set)[g_tables->sample_bv(f)];
//         for (auto idx:idxs) {
//             results.push_back(nameProb[idx].first);
//         }
//     }
// }
// int mdSubsetSampling::add(WeightType addWeight, CountType idx) {}
// int mdSubsetSampling::del(CountType delIdx) {}