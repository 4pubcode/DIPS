#pragma once
#include <string.h>

#include <chrono>
#include <random>

#include "robin_hood.h"
#include "symbols.h"
// #include "../SFMT/dSFMT/dSFMT.h"
// #include "../SFMT/dSFMT/dSFMT.c"

using namespace std;
// std::random_device rd;
// std::default_random_engine generator(rd());
// std::minstd_rand generator(rd());
// std::uniform_real_distribution<double> advanced_unif(0.0, 1.0);
#define BASE 4

double mylog(double x) {
    if (BASE == 4) {
        return std::log(x) / std::log(4);
    } else if (BASE == 3) {
        return std::log(x) / std::log(3);
    } else if (BASE == 2) {
        return std::log2(x);
    } else {
        exit(-1);
    }
}

double mylog(size_t x) {
    if (BASE == 4) {
        return std::log(x) / std::log(4);
    } else if (BASE == 3) {
        return std::log(x) / std::log(3);
    } else if (BASE == 2) {
        return std::log2(x);
    } else {
        exit(-1);
    }
}

double mylog(int x) {
    if (BASE == 4) {
        return std::log(x) / std::log(4);
    } else if (BASE == 3) {
        return std::log(x) / std::log(3);
    } else if (BASE == 2) {
        return std::log2(x);
    } else {
        exit(-1);
    }
}

double mylog(long long x) {
    if (BASE == 4) {
        return std::log(x) / std::log(4);
    } else if (BASE == 3) {
        return std::log(x) / std::log(3);
    } else if (BASE == 2) {
        return std::log2(x);
    } else {
        exit(-1);
    }
}

long long mypow(int n) {
    if (BASE == 4) {
        return (long long)1 << (2 * n);
    } else if (BASE == 3) {
        return pow(3, n);
    } else if (BASE == 2) {
        return pow(2, n);
    } else {
        exit(-1);
    }
}

double unif() {
    return dsfmt_gv_genrand_close_open();
    // return advanced_unif(generator);
    // return ((double)rand() / (RAND_MAX));
}

double geo(double p) {
    // std::geometric_distribution<int> advanced_geo(p);
    // return advanced_geo(generator);
    // return floor(log((double)rand() / (RAND_MAX))/log(1-p));
    return floor(log(unif()) / log(1 - p));
}

double geo(double p, double q) {
    return floor(log(1-q*unif()) / log(1 - p));
}

long long power(long long m, long long n) {
    if (n == 0)
        return 1;
    else if (n % 2 == 0)
        return power(m, n / 2) * power(m, n / 2);
    else
        return m * power(m, n / 2) * power(m, n / 2);
}

int getRi(double weight) {
    // get range index
    // return static_cast<int>(std::ceil(mylog(weight)));
    return static_cast<int>(mylog(weight));
}

int getRi(long long weight) {
    // get range index
    // return static_cast<int>(std::ceil(mylog(weight)));
    return static_cast<int>(mylog(weight));
}

int getIi(int ri, int logn) {
    // get interval index
    return ri / logn;
}

int encode(vector<WeightType> &vec, int A) {
    int m = vec.size();
    int f = 0;
    for (int i = m - 1; i >= 0; i--) {
        f = f * A + static_cast<int>(vec[i]);
    }
    return f;
}

int encode(vector<int> &vec, int A) {
    int m = vec.size();
    int f = 0;
    for (int i = m - 1; i >= 0; i--) {
        f = f * A + vec[i];
    }
    return f;
}

// class Timer {
//  public:
//   Timer(const std::string &task) : _task(task) {
//     std::cout << "-----" << _task << " starts!" << std::endl;
//     _start = std::chrono::system_clock::now();
//   };

//   ~Timer() {
//     _end = std::chrono::system_clock::now();
//     std::cout << "-----" << _task << " ends! " << "Time elapsed: "
//               << std::chrono::duration_cast<std::chrono::microseconds>(_end -
//                                                                        _start)
//                      .count()
//               << " [us]"
//               << std::endl;
//   }

//  private:
//   std::string _task;
//   std::chrono::system_clock::time_point _start;
//   std::chrono::system_clock::time_point _end;
// };

int base_m_upd(int src, int i, int v, int m) {
    // Calculate the value of the digit at position i
    if (v < 0) {
        std::cout << v << std::endl;
        exit(-1);
    }
    int current_digit_value = (src / static_cast<int>(pow(m, i))) % m;

    // Calculate the number with the i-th digit set to 0
    int src_zeroed = src - current_digit_value * static_cast<int>(pow(m, i));

    // Calculate the new number by adding the new digit value at position i
    int new_value = src_zeroed + v * static_cast<int>(pow(m, i));

    return new_value;
}

// int base_m_get(int src, int i, int m) {
//   return (src / static_cast<int>(pow(m, i))) % m;;
// }
