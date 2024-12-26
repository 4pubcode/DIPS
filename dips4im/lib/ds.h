#pragma once
#include <cmath>
#include <iostream>
#include <numeric>
#include <thread>
#include <unordered_map>
#include <vector>

#include "lookup_table.h"
#include "range.h"
#include "robin_hood.h"
// #include "symbols.h"
#include "util.h"

#define CACHE_IN_DOUBLES 16
// #define RI_MAX 512
// #define II_MAX 52
#define MAX_RI 32
using namespace std;

extern lookup_table *g_tables;
extern vector<vector<int>> *g_bits2set;

typedef struct my_index my_index;
typedef struct interval interval;

struct interval {
   private:
    int Ii;

   public:
    vector<WeightType> nweights;
    robin_hood::unordered_map<int, int> Ri2ind;
    // vector<WeightType> aggressive_arr;
    // int r;
    vector<CountType> ris;
    my_index *m_idx;
    int level;
    int f;
    WeightType W;
    WeightType nW;
    WeightType nWf;
    int offset;
    int logn;

   public:
    interval()
        : Ii(0),
          nweights(vector<WeightType>()),
          Ri2ind(robin_hood::unordered_map<int, int>()),
          ris(vector<CountType>()),
          m_idx(nullptr),
          level(0),
          f(0),
          W(0),
          nW(0),
          nWf(0),
          offset(0),
          logn(0) {}
    interval(int Ii, int offset, int logn, int level, int f)
        : Ri2ind(robin_hood::unordered_map<int, int>()),
          // Ri2ind(vector<int>(MAX_RI,-1)),
          //   aggressive_arr(vector<WeightType>(16)),
          //   r(-1),
          ris(vector<CountType>()),
          nweights(vector<WeightType>()),
          m_idx(nullptr),
          Ii(Ii),
          offset(offset),
          logn(logn),
          level(level),
          f(f),
          W(0),
          nW(0),
          nWf(0) {}

    ~interval() {
        if (m_idx) {
            delete m_idx;
        }
        Ri2ind.clear();
        ris.clear();
        nweights.clear();
    }

    inline size_t size(void) const { return ris.size(); }

    void ins(int Ri, WeightType w);

    void del(int Ri, WeightType w);

    void build_my_index();

    void naive_sample(vector<CountType> &ret) const {
        // TODO: direct access results, reduce returns
        //  vector<int> ret;
        //  for (int i = 0; i < arr.size(); i++) {
        //      if (dsfmt_gv_genrand_open_close() * nW <= nweights[i]) {
        //          ret.push_back(arr[i]);
        //      }
        //  }
        //  return ret;
        int n = ris.size();

        // int block=4;
        // for (int i=0; i<ris.size()/block; i++) {
        //     auto w
        // }
        for (int i = 0; i < ris.size(); i++) {
            if (dsfmt_gv_genrand_open_close() * nW <= nweights[i]) {
                ret.push_back(ris[i]);
            }
        }
        // }
        return;
    }
};

struct my_index {
   private:
    robin_hood::unordered_map<CountType, int> e2ind;
    // robin_hood::unordered_map<int, int> e2Ri;

    int init_n;
    // int f;
    int logn;
    int level;

   public:
    vector<CountType> elements;
    vector<WeightType> weights;
    int bv;
    WeightType W;
    robin_hood::unordered_map<int, interval>
        Ii2itv;                                  // Interval ind to interval
    robin_hood::unordered_map<int, range> Ri2r;  // Range ind to range

   public:
    my_index(vector<CountType> &elements, vector<WeightType> &weights, int logn,
             int level);
    ~my_index() {}
    void ins(CountType e, WeightType w);
    void _del(CountType e, WeightType w);
    void del(CountType e);

    void mysample(vector<CountType> &results);

    void sample(vector<CountType> &results);
    void sample_with_interval(const interval &itv, vector<CountType> &results);
    void print_status(void);
};

void my_index::sample_with_interval(const interval &itv,
                                    vector<CountType> &results) {
    if (itv.size() < CACHE_IN_DOUBLES) {
        int n = itv.ris.size();
        for (int i = 0; i < itv.ris.size(); i++) {
            if (dsfmt_gv_genrand_open_close() * itv.nW <= itv.nweights[i]) {
                Ri2r[itv.ris[i]].sample(results, level == 1);
            }
        }
    } else {
        if (level == 1) {
            vector<CountType> index_results;
            itv.m_idx->sample(index_results);
            for (const auto &idx : index_results) {
                Ri2r[idx].sample(results);
            }
        } else {
            auto idxs = (*g_bits2set)[g_tables->sample_bv(itv.f)];
            for (auto idx:idxs) {
                Ri2r[idx].sample(results);
            }
        }
    }
    return;
}

void my_index::mysample(vector<CountType> &results) {
    results.resize(0);
    if (elements.size() == 1) {
        results.push_back(elements[0]);
        return;
    }
    if (elements.size() <= 8) {
        for (int i = 0; i < elements.size(); i++) {
            if (dsfmt_gv_genrand_open_close() * W <= MY_C * weights[i]) {
                results.push_back(elements[i]);
            }
        }
        return;
    }
    if (Ri2r.size() == 1) {
        Ri2r.begin()->second.sample(results, true);
        return;
    }
    if (Ii2itv.size() == 1) {
        const auto &interval = Ii2itv.begin()->second;
        if (interval.m_idx == nullptr) {
            Ri2r[interval.ris[0]].sample(results, true);
        } else {
            sample_with_interval(interval, results);
        }
        return;
    }

    int n = elements.size();
    int log_n = static_cast<int>(ceil(mylog(n)));
    int r = static_cast<int>(static_cast<int>(log2(bv)) / ceil(mylog(n)));
    WeightType remaining_weight = W;

    results.reserve(BASE);  // Reserve space based on the size of elements

    for (int i = r; i >= r - 2 && Ii2itv.find(i) != Ii2itv.end(); --i) {
        const auto &interval = Ii2itv[i];
        remaining_weight -= interval.W;
        // if (dsfmt_gv_genrand_open_close() * W <= interval.W) {
        //     if (interval.m_idx == nullptr) {
        //         Ri2r[interval.ris[0]].sample(results, true);
        //     } else {
        //         sample_with_interval(interval, results);
        //     }
        // }
        auto temp_results = vector<CountType>();
        if (interval.m_idx == nullptr) {
            Ri2r[interval.ris[0]].sample(temp_results, true);
        } else {
            sample_with_interval(interval, temp_results);
        }
        for (auto e: temp_results) {
            if (dsfmt_gv_genrand_open_close() * W <= interval.W) {
                results.push_back(e);
            }
        }  
    }

    if (dsfmt_gv_genrand_open_close() * W <= remaining_weight) {
        for (size_t i = 0; i < elements.size(); ++i) {
            if (weights[i] < mypow((r - 2) * log_n-1) &&
                dsfmt_gv_genrand_open_close() * remaining_weight <=
                    MY_C * weights[i]) {
                results.push_back(elements[i]);
            }
        }
    }

    return;
}

void my_index::print_status(void) {
    for (auto iw : Ii2itv) {
        cout << iw.first << " " << iw.second.W << endl;
    }
    cout << "-----------------" << endl;
}

void my_index::sample(vector<CountType> &results) {
    if (Ii2itv.size() == 1) {
        const auto &interval = Ii2itv[0];
        if (interval.m_idx == nullptr) {
            Ri2r[interval.ris[0]].sample(results);
        } else {
            // vector<CountType> index_results;
            // if (interval.size() < CACHE_IN_DOUBLES) {
            //     interval.naive_sample(index_results);
            // } else {
            //     interval.m_idx->sample(index_results);
            // }
            // // const auto &index_results = interval.m_idx->sample();
            // for (const auto &idx : index_results) {
            //     Ri2r[idx].sample(results);
            // }
            sample_with_interval(interval, results);
        }
        return;
    }
    int r = static_cast<int>(static_cast<int>(log2(bv)) / logn);
    WeightType remaining_weight = W;

    results.reserve(BASE);  // Reserve space based on the size of elements

    for (int i = r; i >= r - 2 && Ii2itv.find(i) != Ii2itv.end(); --i) {
        const auto &itv = Ii2itv.at(i);
        remaining_weight -= itv.W;
        // if (dsfmt_gv_genrand_open_close() * W <= itv.W) {
        //     if (itv.size() == 1) {
        //         Ri2r.at(itv.ris[0]).sample(results);
        //     } else {
        //         sample_with_interval(itv, results);
        //     }
        // }
        auto temp_results = vector<CountType>();
        if (itv.size() == 1) {
            Ri2r.at(itv.ris[0]).sample(temp_results);
        } else {
            sample_with_interval(itv, temp_results);
        }
        for (auto e: temp_results) {
            if (dsfmt_gv_genrand_open_close() * W <= itv.W) {
                results.push_back(e);
            }
        }
    }

    if (dsfmt_gv_genrand_open_close() * W <= remaining_weight) {
        for (size_t i = 0; i < elements.size(); ++i) {
            if (weights[i] <= mypow((r - 3) * logn) &&
                dsfmt_gv_genrand_open_close() * remaining_weight <=
                    weights[i]) {
                results.push_back(elements[i]);
            }
        }
    }

    return;
}

my_index::my_index(vector<CountType> &elements, vector<WeightType> &weights,
                   int logn, int level)
    : bv(0),
      elements(elements),
      weights(weights),
      e2ind(robin_hood::unordered_map<CountType, int>()),
      //   e2Ri(robin_hood::unordered_map<int, int>()),
      Ri2r(robin_hood::unordered_map<int, range>()),
      Ii2itv(robin_hood::unordered_map<int, interval>()),
      level(level),
      logn(logn) {
    W = 0;
    int cnt = 0;
    for (int i = 0; i < weights.size(); i++) {
        auto e = elements[i];
        auto w = weights[i];
        cnt += sizeof(CountType) + sizeof(WeightType);
        auto Ri = getRi(w);
        auto Ii = getIi(Ri, logn);
        // if (f == -1) {
        e2ind[e] = i;
        // e2Ri[e] = Ri;
        // }
        W += w;
        if (Ri2r.find(Ri) == Ri2r.end()) {
            bv |= (1 << Ri);
            Ri2r[Ri] = range(mypow(Ri + 1));
            Ri2r[Ri].ins(e, w);
            cnt += sizeof(CountType) + sizeof(WeightType);
            if (Ii2itv.find(Ii) == Ii2itv.end()) {
                // Ii2itv[Ii] = interval(Ii, Ii * logn);
                Ii2itv[Ii] =
                    interval(Ii, Ii * logn, logn, level, (level == 2) ? 0 : -1);
                // TODO:use w or Ri2r[Ri].W?
                Ii2itv[Ii].ins(Ri, w);
            } else {
                Ii2itv[Ii].ins(Ri, w);
            }
        } else {
            Ri2r[Ri].ins(e, w);
            // TODO:is this safe? will there always be a I for this R in this
            // case?
            Ii2itv[Ii].ins(Ri, w);
        }
    }
}

void my_index::ins(CountType e, WeightType w) {
    // checked: 11 Jul 5:36 PM
    // elements, weights, e2ind
    // init_n, f, logn
    // bv, W, Ii2itv, Ri2r
    int ind_e;
    W += w;
    if (e2ind.find(e) == e2ind.end()) {  // e was not in the problem instance
        e2ind[e] = weights.size();
        ind_e = weights.size();
        weights.push_back(w);
        elements.push_back(e);

        auto Ri = getRi(w);
        auto Ii = getIi(Ri, logn);

        if (Ri2r.find(Ri) == Ri2r.end()) {  // the range has no other elements
            Ri2r[Ri] = range(mypow(Ri + 1));
            bv |= (1 << Ri);
            Ri2r[Ri].ins(e, w);
            if (Ii2itv.find(Ii) ==
                Ii2itv.end()) {  // the interval has no other ranges
                Ii2itv[Ii] =
                    interval(Ii, Ii * logn, logn, level, (level == 2) ? 0 : -1);
            }
            Ii2itv[Ii].ins(Ri, w);
        } else {  // the range has other elements
            Ri2r[Ri].ins(e, w);
            Ii2itv[Ii].ins(Ri, w);
        }
    } else {  // e is currently in the problem instance
        ind_e = e2ind[e];
        // snapshot before update
        auto old_we = weights[ind_e];
        auto oldRi = getRi(weights[ind_e]);
        // update
        weights[ind_e] += w;
        auto Ri = getRi(weights[ind_e]);
        auto Ii = getIi(Ri, logn);
        if (oldRi == Ri) {  // e is in the same range as before
            Ri2r[Ri].ins(e, w);
            Ii2itv[Ii].ins(Ri, w);
        } else {  // e should be put in another range
            Ri2r[oldRi].del(e, old_we);
            if (Ri2r[oldRi].size() == 0) {
                Ri2r.erase(oldRi);
                bv ^= (1 << oldRi);
            }

            auto oldIi = getIi(oldRi, this->logn);
            Ii2itv[oldIi].del(oldRi, old_we);
            if (Ii2itv[oldIi].size() == 0) {
                Ii2itv.erase(oldIi);
            }

            if (Ri2r.find(Ri) ==
                Ri2r.end()) {  // the range has no other elements
                Ri2r[Ri] = range(mypow(Ri + 1));
                bv |= (1 << Ri);
                Ri2r[Ri].ins(e, weights[ind_e]);
                if (Ii2itv.find(Ii) ==
                    Ii2itv.end()) {  // the interval has no other ranges
                    Ii2itv[Ii] = interval(Ii, Ii * logn, logn, level,
                                          (level == 2) ? 0 : -1);
                }
                Ii2itv[Ii].ins(Ri, weights[ind_e]);
            } else {  // the range has other elements
                Ri2r[Ri].ins(e, weights[ind_e]);
                Ii2itv[Ii].ins(Ri, weights[ind_e]);
            }

            // Ri2r[Ri].ins(e, weights[ind_e]);
            // auto oldIi = getIi(oldRi, this->logn);
            // Ii2itv[oldIi].del(oldRi, old_we);
            // if (Ii2itv[oldIi].W == 0) {
            //     Ii2itv.erase(oldIi);
            // }
            // Ii2itv[Ii].ins(Ri, weights[ind_e]);
        }
    }
}

void my_index::_del(CountType e, WeightType w) {
    // checked (not confidently): 11 Jul 6:01 PM
    // elements, weights, e2ind
    // init_n, f, logn
    // bv, W, Ii2itv, Ri2r
    auto ind_e = e2ind[e];
    W -= w;
    if (W < 0) {
        auto a = 0;
    }

    // snapshot before update
    auto old_we = weights[ind_e];
    auto oldRi = getRi(weights[ind_e]);
    auto oldIi = getIi(oldRi, logn);
    // update
    weights[ind_e] -= w;
    if (weights[ind_e] <= 0) {  // e should be cleared
        elements[ind_e] = elements[elements.size() - 1];
        elements.pop_back();
        e2ind[elements[ind_e]] = ind_e;
        e2ind.erase(e);

        weights[ind_e] = weights[weights.size() - 1];
        weights.pop_back();

        Ri2r[oldRi].del(e, w);
        if (Ri2r[oldRi].size() == 0) {
            Ri2r.erase(oldRi);
            bv ^= (1 << oldRi);
        }
        Ii2itv[oldIi].del(oldRi, w);
    } else {  // e still exists
        auto Ri = getRi(weights[ind_e]);
        auto Ii = getIi(Ri, logn);
        if (oldRi == Ri) {  // e is in the same range as before
            Ri2r[Ri].del(e, w);
            Ii2itv[Ii].del(Ri, w);
        } else {  // e should be put in another range
            Ri2r[oldRi].del(e, old_we);
            if (Ri2r[oldRi].size() == 0) {
                Ri2r.erase(oldRi);
                bv ^= (1 << oldRi);
            }

            Ii2itv[oldIi].del(oldRi, old_we);
            if (Ii2itv[oldIi].size() == 0) {
                Ii2itv.erase(oldIi);
            }

            if (Ri2r.find(Ri) ==
                Ri2r.end()) {  // the range has no other elements
                Ri2r[Ri] = range(mypow(Ri + 1));
                bv |= (1 << Ri);
                Ri2r[Ri].ins(e, weights[ind_e]);
                if (Ii2itv.find(Ii) ==
                    Ii2itv.end()) {  // the interval has no other ranges
                    Ii2itv[Ii] = interval(Ii, Ii * logn, logn, level,
                                          (level == 2) ? 0 : -1);
                }
                Ii2itv[Ii].ins(Ri, weights[ind_e]);
            } else {  // the range has other elements
                Ri2r[Ri].ins(e, weights[ind_e]);
                Ii2itv[Ii].ins(Ri, weights[ind_e]);
            }

            // Ri2r[Ri].ins(e, weights[ind_e]);
            // Ii2itv[oldIi].del(oldRi, old_we);
            // Ii2itv[Ii].ins(Ri, weights[ind_e]);
        }
    }
    if (Ii2itv[oldIi].size() == 0) {
        Ii2itv.erase(oldIi);
    }
}

void my_index::del(CountType e) { _del(e, weights[e2ind[e]]); }

void interval::build_my_index() {
    // remember to normalize
    // vector<WeightType> *nweights = new vector<WeightType>(weights.size());
    // for (size_t i = 0; i < weights.size(); i++) {
    //     // TODO:fix WeightType /WeightType
    //     (*nweights)[i] = weights[i] / (WeightType)mypow(offset);
    // }
    m_idx = new my_index(ris, nweights, static_cast<int>(ceil(mylog(logn))),
                         level + 1);
    auto a = 0;
}

void interval::ins(int Ri, WeightType w) {
    auto nw = w / (WeightType)mypow(offset);
    nW += nw;
    W += w;

    if (Ri2ind.find(Ri) == Ri2ind.end()) {  // Ri hasn't been created
        Ri2ind[Ri] = ris.size();
        ris.push_back(Ri);
        nweights.push_back(nw);
        if ((m_idx == nullptr) && (Ri2ind.size() > 1) && (f == -1)) {
            build_my_index();
        }
        nWf += floor(nw);
    } else {  // Ri was created before
        auto ind_ri = Ri2ind[Ri];
        auto old_nw = nweights[ind_ri];
        nWf -= floor(old_nw);
        nweights[ind_ri] += nw;
        nWf += floor(nweights[ind_ri]);
        if (m_idx != nullptr) {  // update index
            // if (getRi(old_nw) == getRi(nweights[ind_ri])) {
            //     m_idx->ins(Ri, nw);
            // } else { // the total weight of Ri change greatly
            //     m_idx->del(Ri, old_nw);
            //     m_idx->ins(Ri, nweights[ind_ri]);
            // }
            m_idx->ins(Ri, nw);
        }
    }
    if (f != -1) {
        // f = base_m_upd(f, Ri - offset, weights[Ri2loc[Ri]] /
        // mypow(offset),
        //                pow(mypow(logn), 2));
        f = base_m_upd(f, Ri - offset, ceil(nweights[Ri2ind[Ri]]),
                       pow(mypow(logn), 2));
    }
}

void interval::del(int Ri, WeightType w) {
    // nweights, Ri2ind, ris, m_idx
    // f, W, nW, nWf, offset, logn
    // Ii
    auto ind_ri = Ri2ind[Ri];
    auto nw = w / mypow(offset);
    W -= w;
    if (W < 0) {
        auto a = 0;
    }
    nW -= nw;

    nWf -= floor(nweights[ind_ri]);
    auto old_nw = nweights[ind_ri];
    nweights[ind_ri] -= nw;
    nWf += floor(nweights[ind_ri]);

    if (f != -1) {
        // f = base_m_upd(f, Ri - offset, weights[Ri2loc[Ri]] /
        // mypow(offset),
        //                pow(mypow(logn), 2));
        f = base_m_upd(f, Ri - offset, ceil(nweights[ind_ri]),
                       pow(mypow(logn), 2));
    }

    if (nweights[ind_ri] <= 0) {  // considering precision
        ris[ind_ri] = ris[ris.size() - 1];
        ris.pop_back();
        Ri2ind[ris[ind_ri]] = ind_ri;
        Ri2ind.erase(Ri);

        nweights[ind_ri] = nweights[nweights.size() - 1];
        nweights.pop_back();
        if (Ri2ind.size() == 1) {
            delete m_idx;
            m_idx = nullptr;
        }
    }

    if (m_idx != nullptr) {  // update index
        // if (getRi(old_nw) == getRi(nweights[ind])) {
        //     m_idx->del(Ri, nw);
        //     m_idx->ins(Ri, nweights[ind]);
        // } else {
        //     m_idx->del(Ri, old_nw);
        //     m_idx->ins(Ri, nweights[ind]);
        // }
        m_idx->_del(Ri, nw);
        if (m_idx->W == 0) {
            exit(-1);
        }
    }
}
