#pragma once
#include <algorithm>
#include <cstring>
#include <stack>

#include "symbols.h"
using namespace std;

class Alias {
   public:
    double *p;
    CountType *h;
    CountType *map1;
    CountType n;
    Alias() {}
    void input(pair<CountType, double> *pi, CountType size) {
        n = size;
        double sum = 0;
        // stack<CountType> small;
        // stack<CountType> big;
        // 初始化代价比stack高
        CountType *small = new CountType[n];
        CountType *big = new CountType[n];
        CountType small_cnt = 0, big_cnt = 0;
        p = new double[n];
        h = new CountType[n];
        map1 = new CountType[n];
        for (CountType i = 0; i < n; i++) {
            sum += pi[i].second;
            map1[i] = pi[i].first;
        }
        for (CountType i = 0; i < n; i++) {
            p[i] = pi[i].second * n / sum;
            if (p[i] > 1)
                // big.push(i);
                big[big_cnt++] = i;
            else
                // small.push(i);
                small[small_cnt++] = i;
        }
        // cnt的意义变成指示
        small_cnt--;
        big_cnt--;
        while (small_cnt >= 0 && big_cnt >= 0 && small_cnt < n && big_cnt < n) {
            CountType smallId = small[small_cnt];
            CountType bigId = big[big_cnt];
            h[smallId] = bigId;
            p[bigId] -= (1 - p[smallId]);
            if (p[bigId] < 1) {
                small[small_cnt] = bigId;
                big_cnt--;
            } else
                small_cnt--;
        }
        delete[] small;
        delete[] big;
    }

    Alias(pair<CountType, double> *pi, CountType size) {
        n = size;
        double sum = 0;
        // stack<CountType> small;
        // stack<CountType> big;
        // 初始化代价比stack高
        CountType *small = new CountType[n];
        CountType *big = new CountType[n];
        CountType small_cnt = 0, big_cnt = 0;
        p = new double[n];
        h = new CountType[n];
        map1 = new CountType[n];
        for (CountType i = 0; i < n; i++) {
            sum += pi[i].second;
            map1[i] = pi[i].first;
        }
        for (CountType i = 0; i < n; i++) {
            p[i] = pi[i].second * n / sum;
            if (p[i] > 1)
                // big.push(i);
                big[big_cnt++] = i;
            else
                // small.push(i);
                small[small_cnt++] = i;
        }
        // cnt的意义变成指示
        small_cnt--;
        big_cnt--;
        while (small_cnt >= 0 && big_cnt >= 0) {
            CountType smallId = small[small_cnt];
            CountType bigId = big[big_cnt];
            h[smallId] = bigId;
            p[bigId] -= (1 - p[smallId]);
            if (p[bigId] < 1) {
                small[small_cnt] = bigId;
                big_cnt--;
            } else
                small_cnt--;
        }
        delete[] small;
        delete[] big;
    }

    Alias(const Alias &tmp) {
        n = tmp.n;
        p = new double[n];
        h = new CountType[n];
        map1 = new CountType[n];
        memcpy(p, tmp.p, sizeof(double) * n);
        memcpy(h, tmp.h, sizeof(CountType) * n);
        memcpy(map1, tmp.map1, sizeof(CountType) * n);
    }

    Alias &operator=(const Alias &tmp) {
        if (&tmp != this) {
            n = tmp.n;
            p = new double[n];
            h = new CountType[n];
            map1 = new CountType[n];
            memcpy(p, tmp.p, sizeof(double) * n);
            memcpy(h, tmp.h, sizeof(CountType) * n);
            memcpy(map1, tmp.map1, sizeof(CountType) * n);
        }

        return *this;
    }

    ~Alias() {
        delete[] p;
        delete[] h;
        delete[] map1;
    }

    CountType generateRandom(const double &rand1, const double &rand2) {
        CountType firstId = rand1 * n;
        CountType answer =
            rand2 < p[firstId] ? map1[firstId] : map1[h[firstId]];
        return answer;
    }
};
