#include "lookup_table.h"
#include "ds.h"
#include "range.h"
#include "ppsMethod.h"
// #define BASE 4
#define M 3.0
// #define RI_MAX 100
// #define II_MAX 10

extern lookup_table *g_tables;
extern vector<vector<int>> *g_bits2set;

class DIPS : public ppsMethod {
   private:
    my_index *level1;

    // lookup_table *tables;
    // vector<vector<int>> *bits2set;
    // int n;
    int oldSize;
    WeightType W;
    int f;
    vector<CountType> es;
    vector<WeightType> ws;
    unordered_map<CountType, int> e2ind;

   public:
    DIPS(/* args */);
    ~DIPS();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &results);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
};

DIPS::DIPS(/* args */) : level1(nullptr), W(0), f(0) {}

DIPS::~DIPS() {
    if (level1) {
        delete level1;
    }
    // delete level1;
}

int DIPS::init(vector<CountType> elements,
                           vector<WeightType> weights) {
    // n = elements.size();
    oldSize = elements.size();
    int logn = static_cast<int>(ceil(mylog(elements.size())));
    // TODO this is bottleneck, and space is too large
    // es = elements;
    // ws = weights;
    level1 = new my_index(elements, weights, max(1, logn), 1);
    // level1->print_status();
    for (auto &couple : level1->Ii2itv) {
        if (couple.second.size() != 1) {
            couple.second.build_my_index();
            // couple.second.m_idx->print_status();
        } else {
            couple.second.m_idx = nullptr;
        }
    }
    
    return 1;
}

int DIPS::query(vector<CountType> &results) {
    level1->mysample(results);
    return 1;
}

int DIPS::add(WeightType addWeight, CountType idx) {
    level1->ins(idx, addWeight);
    if (level1->elements.size() >=2*oldSize) {
        init(level1->elements, level1->weights);
    }
    return 1;
}

int DIPS::del(CountType delIdx) {
    level1->del(delIdx);
    if (level1->elements.size() <= oldSize/2) {
        init(level1->elements, level1->weights);
    }
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

// class DIPS : public subsetSamplingMethod {
//    private:
//     CountType n;
//     CountType size;
//     pair<CountType, WeightType> *nameProb;
//     CountType *dict;
//     vector<CountType> *ranges;
//     WeightType *wr;
//     DIPS **intervals;
//     WeightType *wi;

//     // lookup_table *tables;
//     // vector<vector<int>> *bits2set;
//     WeightType W;
//     int f;

//    public:
//     DIPS(/* args */);
//     ~DIPS();
//     int init(vector<CountType> elements, vector<WeightType> weights);
//     int query(vector<CountType> &results);
//     int add(WeightType addWeight, CountType idx);
//     int del(CountType delIdx);
// };

// int DIPS::init(vector<CountType> elements,
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
//     intervals = new DIPS*[II_MAX];
//     wi = new WeightType[II_MAX];
//     for (CountType i=0; i<II_MAX; i++) {
//         intervals[i] = nullptr;
//         wi[i] = 0;
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
// int DIPS::query(vector<CountType> &results) {
//     if (n < M) {
//         auto idxs = (*g_bits2set)[g_tables->sample_bv(f)];
//         for (auto idx:idxs) {
//             results.push_back(nameProb[idx].first);
//         }
//     }
// }
// int DIPS::add(WeightType addWeight, CountType idx) {}
// int DIPS::del(CountType delIdx) {}