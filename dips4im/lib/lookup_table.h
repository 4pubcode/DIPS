#pragma once
#include "util.h"

class lookup_table {
    // m <= 5 in practice
   private:
    vector<vector<uint8_t>> tables;
    int m;
    int base;
    int A;

   public:
    lookup_table(int m, int base) : m(m), base(base), A(pow(mypow(m), 2)) {
        tables = vector<vector<uint8_t>>(pow(A, m), vector<uint8_t>());
        // tables = vector<vector<uint8_t>>(pow(A, m));
        buildLookupTables(m, A);
    }

    inline uint8_t sample_bv(int f) {
        // return tables[f][rand() % tables[f].size()];
        return tables[f][static_cast<int>(unif() * tables[f].size())];
    }

   private:
    void buildLookupTable(vector<uint8_t> &table, const vector<int> &weights);

    void buildLookupTables(int m, int A);

    void findCombinations(int n, int m, int nnz, int start,
                          vector<int> &current, vector<vector<int>> &result);

    vector<vector<int>> sumCombinations(int n, int m);

    void generateSums(const vector<int> &partition, int index,
                      vector<int> &current);

    void generatePartitionSums(const vector<int> &partition);
};

void lookup_table::buildLookupTable(vector<uint8_t> &table,
                                    const vector<int> &weights) {
    int W = accumulate(weights.begin(), weights.end(), 0) - m;
    int totalEvents = 1 << m;  // 2^m events

    for (int event = 0; event < totalEvents; ++event) {
        int eventWeight = 1;
        uint8_t eventOutcome = 0;

        for (int i = 0; i < m; ++i) {
            if (event & (1 << i)) {  // if the i-th bit is set
                eventWeight *= weights[i];
                eventOutcome |= (1 << i);  // little-endian
            } else {
                eventWeight *= (W - weights[i]);
                // eventOutcome = eventOutcome << 1;
            }
        }

        // table.insert(table.begin(), eventWeight, eventOutcome);
        for (int j = 0; j < eventWeight; j++) {
            table.push_back(eventOutcome);
        }
    }
}

void lookup_table::buildLookupTables(int m, int A) {
    for (int i = 1; i <= mypow(m); i++) {
        vector<vector<int>> combinations = sumCombinations(i, m);
        for (const auto &partition : combinations) {
            generatePartitionSums(partition);
        }
    }
}

// Helper function to find all combinations
void lookup_table::findCombinations(int n, int m, int nnz, int start,
                                    vector<int> &current,
                                    vector<vector<int>> &result) {
    if (m == 1) {
        current.push_back(n);
        if ((nnz + (n != 0)) > 1) {
            result.push_back(current);
        }
        current.pop_back();
        return;
    }

    for (int i = start; i <= n; ++i) {
        current.push_back(i);
        findCombinations(n - i, m - 1, nnz + (i != 0), 0, current, result);
        current.pop_back();
    }
}

// Main function to find all possible ways to represent n as the sum of m
// nonnegative integers
vector<vector<int>> lookup_table::sumCombinations(int n, int m) {
    vector<vector<int>> result;
    vector<int> current;
    findCombinations(n, m, 0, 0, current, result);
    return result;
}

// Helper function to recursively generate the possible sums
void lookup_table::generateSums(const vector<int> &partition, int index,
                                vector<int> &current) {
    if (index == partition.size()) {
        // result.push_back(current);
        auto f = encode(current, A);
        if (tables[f].empty()) {
            buildLookupTable(tables[f], current);
        }
        return;
    }

    int range_start = partition[index] * mypow(index);
    int range_end = partition[index] * mypow(index + 1);

    for (int j = range_start; j <= range_end; ++j) {
        current[index] = j;
        generateSums(partition, index + 1, current);
    }
}

// Main function to generate all possible sums based on the partition and
// base
void lookup_table::generatePartitionSums(const vector<int> &partition) {
    vector<int> current(partition.size(), 0);
    generateSums(partition, 0, current);
}