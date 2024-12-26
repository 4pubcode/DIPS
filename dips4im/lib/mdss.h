#pragma once
#include <math.h>

#include <algorithm>
#include "robin_hood.h"
#include <unordered_map>
#include <vector>
#include <utility>
#include <tuple>
#include "util.h"
using namespace std;
typedef unsigned int uint;
typedef uint CountType;
typedef double WeightType;
class md_index {
 private:
  vector<vector<pair<CountType,WeightType>>> buckets1;
  WeightType rem1q;
  vector<CountType> rem1;
  vector<vector<CountType>> buckets2;
  vector<CountType> rem2;
  WeightType rem2q;
  vector<WeightType> e2p2;
  vector<CountType> e2;
  vector<CountType> e3;
  vector<WeightType> e2p3;
  vector<CountType> elems;
  vector<WeightType> ps;
  // lookup_table *tables;
  // vector<vector<int>> *bits2set;
  
  int n, logn, loglogn;

 public:
 robin_hood::unordered_map<CountType, WeightType> e2p1;
//  unordered_map<CountType, WeightType> e2p1;
//  vector<WeightType> e2p1;
 WeightType W;
  md_index(vector<CountType> elements, vector<WeightType> weights) {
    n = elements.size();
    if (n == 0) return;
    // e2p1.resize(n);
    logn = 2*static_cast<int>(ceil(log2(n)));
    loglogn = logn > 0 ? 2*static_cast<int>(ceil(log2(logn))) : 0;
    buckets1.resize(logn);
    buckets2.resize(loglogn);
    e2p2.resize(logn);
    e3.resize(0);
    e2p3.resize(loglogn);

    W = 0;
    for (int i = 0; i < n; i++) {
      elems.push_back(elements[i]);
      W += weights[i];
    }

    for (int i = 0; i < logn; i++) {
      e2p2[i] = 0;
    }
    rem1q = 0;
    for (int i = 0; i < n; i++) {
      auto e = elements[i];
      auto p = weights[i] / W;
      ps.push_back(p);
      e2p1[e] = p;
      size_t id = -ceil(log2(p));
      if (id < logn) {
        buckets1[id].push_back(pair<CountType,WeightType>(e,p));
        e2p2[id] = 1 - (1 - e2p2[id]) * (1 - p);
      } else {
        rem1.push_back(e);
        rem1q = 1 - (1 - rem1q) * (1 - p);
      }
    }
    for (int i=0; i<logn; i++) {
      if (e2p2[i]>0) {
        e2.push_back(i);
      }
    }

    for (int i = 0; i < loglogn; i++) {
      e2p3[i] = 0;
    }
    rem2q = 0;
    for (int i = 0; i < logn; i++) {
      auto e = i;
      auto p = e2p2[e];
      size_t id = -ceil(log2(p));
      if (id < loglogn) {
        buckets2[id].push_back(e);
        e2p3[id] = 1 - (1 - e2p3[id]) * (1 - p);
      } else {
        rem2.push_back(e);
        rem2q = 1 - (1 - rem2q) * (1 - p);
      }
    }
    for (int i=0; i<loglogn; i++) {
      if (e2p3[i]>0) {
        e3.push_back(i);
      }
    }
  }

  void sample(vector<CountType> &results) {
    // if (n==1) {
    //   results.push_back(elems[0]);
    //   return;
    // }
    // if (n<=16) {
    //   for (int i = 0; i < n; i++) {
    //     if (dsfmt_gv_genrand_open_close() < ps[i]) {
    //       results.push_back(elems[i]);
    //     }
    //   }
    //   return;
    // }
    vector<int> cand2;
    vector<int> cand1;
    if (e2.size()<8) {
      for (auto e:e2) {
        if (dsfmt_gv_genrand_close_open() <= e2p2[e]) {
          cand1.push_back(e);
        }
      }
      goto CAND1FINISHED;
    }
    for (auto i:e3) {
      if (dsfmt_gv_genrand_open_close() < e2p3[i]) {
        cand2.push_back(i);
      }
    }
    
    
    for (auto id : cand2) {
      auto &arr = buckets2[id];
      WeightType q = e2p3[id];
      if (arr.size() == 1) {
        auto e = arr[0];
        if (dsfmt_gv_genrand_close_open() * q <= e2p2[e]) {
          cand1.push_back(e);
        }
        continue;
      }
      WeightType tmp = 1.0 / (1 << id) / q;
      WeightType p = tmp > 1 ? 1 : tmp;
      int loc = geo(p);
      while (loc < arr.size()) {
        auto e = arr[loc];
        if (dsfmt_gv_genrand_close_open() * p * q <= e2p2[e]) {
          cand1.push_back(e);
        }
        loc += geo(p) + 1;
      }
    }
    if (dsfmt_gv_genrand_close_open() <= rem2q) {
      for (auto e : rem2) {
        if (dsfmt_gv_genrand_close_open() * rem2q <= e2p2[e]) {
          cand1.push_back(e);
        }
      }
    }
CAND1FINISHED:
    results.resize(0);
    for (auto id : cand1) {
      auto &arr = buckets1[id];
      WeightType q = e2p2[id];
      if (arr.size() == 1) {
        auto e = arr[0].first;
        // if (dsfmt_gv_genrand_close_open() * q <= e2p1[e]) {
        if (dsfmt_gv_genrand_close_open() * q <= arr[0].second) {
          results.push_back(e);
        }
        continue;
      }
      WeightType tmp = 1.0 / (1 << id) / q;
      WeightType p = tmp > 1 ? 1 : tmp;
      int loc = geo(p);
      while (loc < arr.size()) {
        auto e = arr[loc].first;
        // if (dsfmt_gv_genrand_close_open() * p * q <= e2p1[e]) {
        if (dsfmt_gv_genrand_close_open() * p * q <= arr[0].second) {
          results.push_back(e);
        }
        loc += geo(p) + 1;
      }
    }
    if (dsfmt_gv_genrand_close_open() <= rem1q) {
      for (auto e : rem1) {
        if (dsfmt_gv_genrand_close_open() * rem1q <=  e2p1[e]) {
          results.push_back(e);
        }
      }
    }
  }
};