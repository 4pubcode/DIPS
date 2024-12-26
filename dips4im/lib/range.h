#pragma once
#include "util.h"

typedef struct range {
   private:
    robin_hood::unordered_map<CountType, int> e2ind;

    vector<WeightType> weights;
    WeightType w_bar;

   public:
    vector<CountType> arr;
    WeightType W;

   public:
    range()
        : e2ind(robin_hood::unordered_map<CountType, int>()),
          arr(vector<CountType>()),
          weights(vector<WeightType>()),
          w_bar(0),
          W(0) {}
    range(WeightType w_bar)
        : e2ind(robin_hood::unordered_map<CountType, int>()),
          arr(vector<CountType>()),
          weights(vector<WeightType>()),
          w_bar(w_bar),
          W(0) {}

    ~range() {}

    int size(void) const { return arr.size(); }
    void ins(CountType e, WeightType w) {
        if (e2ind.find(e) == e2ind.end()) {
            e2ind[e] = arr.size();
            arr.push_back(e);
            weights.push_back(w);
        } else {
            weights[e2ind[e]] += w;
        }
        W += w;
    }

    void del(CountType e, WeightType w) {
        auto ind_e = e2ind[e];
        W -= w;
        weights[ind_e] -= w;
        if (weights[ind_e] <= 0) {
            arr[ind_e] = arr[arr.size() - 1];
            e2ind[arr[ind_e]] = ind_e;
            arr.pop_back();
            weights[ind_e] = weights[weights.size() - 1];
            weights.pop_back();
            e2ind.erase(e);
        }
        // auto ind_e = e2loc[e];
        // W -= weights[ind_e];
        // arr[ind_e] = arr[arr.size() - 1];
        // e2loc[arr[ind_e]] = ind_e;
        // arr.pop_back();
        // weights[ind_e] = weights[weights.size() - 1];
        // weights.pop_back();
        // e2loc.erase(e);
    }

    void sample(vector<CountType> &results, bool is_final = false) const {
        if (W < w_bar) {
            for (int i = 0; i < weights.size(); i++) {  // #elements <= BASE
                if (dsfmt_gv_genrand_close_open() * W <=
                    (is_final ? MY_C : 1) * weights[i]) {
                    results.push_back(arr[i]);
                }
            }
            return;
        }
        WeightType p = w_bar / (WeightType)W;
        int loc = geo(p);

        while (loc < arr.size()) {
            if (dsfmt_gv_genrand_close_open() * w_bar <=
                (is_final ? MY_C : 1) * weights[loc]) {
                results.push_back(arr[loc]);
            }
            loc += geo(p) + 1;
        }
        return;
    }
} range;
