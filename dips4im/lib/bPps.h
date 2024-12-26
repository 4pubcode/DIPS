#include "alias.h"
#include "ppsMethod.h"
#include "robin_hood.h"
class bPps : public ppsMethod {
   private:
    CountType n, lastn;
    // vector<CountType> sortedprob;
    CountType blockNum;
    Alias *aliasStruct;
    CountType maxSubsetID;

    // CountType *sortedIdx;
    WeightType *barProb;
    WeightType *blockDist;
    element *sortedprob;
    CountType newSize;
    pair<CountType, WeightType> *nameProb;
    // robin_hood::unordered_map<CountType, CountType> dict;
    WeightType W;

   public:
    bPps(/* args */);
    ~bPps();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &result);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
    // int find(CountType& idx, WeightType w);
};

// int bringmannSubsetSampling::find(CountType& idx, WeightType w) {
//     if (adjList[idx].first == -1) return -1;

// }

bPps::bPps(/* args */) {}

bPps::~bPps() {}

int bPps::query(vector<CountType> &result) {
    if (n == 0) return 0;
    mt19937 generator(time(NULL));
    result.resize(0);
    CountType idx = 0;
    vector<CountType> tmp_results;
    while (idx < blockNum) {
        CountType bucketIdx = aliasStruct[idx].generateRandom(
            dsfmt_gv_genrand_close_open(), dsfmt_gv_genrand_close_open());
        if (bucketIdx == -1) break;
        CountType pivot = bucketIdx == 0 ? 0 : 1 << bucketIdx;
        CountType startMin = pivot;
        CountType endMax = bucketIdx == 0 ? 2 : startMin * 2;
        if (endMax > n) endMax = n;
        WeightType w = barProb[sortedprob[pivot].first];
        const double log_prob = log2(1 - w);
        WeightType selectProb = blockDist[bucketIdx];
        double uni =
            dsfmt_gv_genrand_close_open() + (1 - selectProb) / selectProb;
        startMin += floor(log2(uni * selectProb) / log_prob);
        if (startMin < endMax) {
            if (dsfmt_gv_genrand_close_open() <
                sortedprob[startMin].second / w) {
                tmp_results.push_back(sortedprob[startMin].first);
            }
            startMin++;
        }
        while (startMin < endMax) {
            double uni_prob = dsfmt_gv_genrand_close_open();
            startMin += floor(log2(uni_prob) / log_prob);
            if (startMin >= endMax) break;
            if (dsfmt_gv_genrand_close_open() <
                sortedprob[startMin].second / w) {
                tmp_results.push_back(sortedprob[startMin].first);
            }
            startMin++;
        }
        idx = bucketIdx + 1;
    }
    for (auto ind : tmp_results) {
        result.push_back(nameProb[ind].first);
    }
    return 0;
}

int bPps::init(vector<CountType> elements, vector<WeightType> weights) {
    n = weights.size();
    if (n == 0) return 0;
    // auto newSize = (n > CountTypeMAX / 2) ? CountTypeMAX : n*2;
    auto newSize = n;
    nameProb = new pair<CountType, WeightType>[newSize];
    W = 0;
    for (CountType i = 0; i < n; i++) {
        W += weights[i];
    }
    for (CountType i = 0; i < n; i++) {
        weights[i] /= W;
        nameProb[i] = make_pair(elements[i], weights[i]);
        // dict[elements[i]] = i;
        elements[i] = i;
    }

    vector<CountType> *buckets = new vector<CountType>[INTBITS];
    lastn = n;
    // partition
    maxSubsetID = ceil(log2(n)) > 1 ? ceil(log2(n)) : 1;
    // cout << "original maxSubsetID: " << maxSubsetID << endl;
    if (maxSubsetID >= INTBITS) maxSubsetID = INTBITS - 1;
    // cout << "maxSubsetID: " << maxSubsetID << endl;

    // if (n > CountTypeMAX / 2)
    //     newSize = CountTypeMAX;
    // else
    //     newSize = n * 2;
    newSize = n;

    barProb = new WeightType[newSize];
    for (size_t i = 0; i < n; i++) {
        if (weights[i] < 0) continue;
        CountType subsetID = -(ceil(log2(weights[i])));
        if (subsetID > maxSubsetID) subsetID = maxSubsetID;

        buckets[subsetID].push_back(elements[i]);
        barProb[i] = 1.0 / (1 << subsetID);
    }

    sortedprob = new element[newSize];
    CountType sortedprobsize = 0;
    // get sorted array
    for (CountType i = 0; i <= maxSubsetID; i++) {
        if (buckets[i].size() == 0) continue;
        vector<CountType> &tmpbucket = buckets[i];
        // cout << "# bucket " << i << " start at " << sortedprobsize << endl;
        for (auto &idx : tmpbucket) {
            sortedprob[sortedprobsize++] = make_pair(idx, weights[idx]);
            // sortedIdx[idx] = sortedprobsize - 1;
        }
    }

    // get block probability
    WeightType firstw = barProb[sortedprob[0].first];
    blockDist = new WeightType[INTBITS];
    CountType sizeblockDist = 0;
    blockDist[sizeblockDist++] =
        1 - (1 - firstw) *
                (1 - firstw);  // sortedprob[0] and sortedprob[1] is a group.
    CountType tmpidx = 2;
    while (tmpidx < n) {
        WeightType w = barProb[sortedprob[tmpidx].first];
        CountType s;
        if (CountTypeMAX - tmpidx < tmpidx)
            s = n - tmpidx;
        else
            s = tmpidx * 2 > n ? n - tmpidx : tmpidx;
        blockDist[sizeblockDist++] = 1 - pow(1 - w, s);
        tmpidx *= 2;
    }

    // get skip probability, construct alias
    blockNum = sizeblockDist;
    aliasStruct = new Alias[blockNum];
    for (CountType i = 0; i < blockNum; i++) {
        pair<CountType, double> *pi =
            new pair<CountType, double>[blockNum - i + 1];
        WeightType failProb = 1 - blockDist[i];
        pi[0] = make_pair(i, blockDist[i]);
        CountType idx = i + 1;
        for (CountType j = 1; j < blockNum - i; j++) {
            pi[j] = make_pair(idx, failProb * blockDist[idx]);
            failProb *= 1 - blockDist[idx];
            idx++;
        }
        pi[blockNum - i] = make_pair(-1, failProb);
        aliasStruct[i].input(pi, blockNum - i + 1);
        delete[] pi;
    }
    delete[] buckets;
    return 0;
}

int bPps::add(WeightType addWeight, CountType idx) {
    vector<CountType> elements(n + 1);
    vector<WeightType> weights(n + 1);
    for (int i = 0; i < n; i++) {
        elements[i] = nameProb[i].first;
        weights[i] = nameProb[i].second * W;
    }
    elements[n] = idx;
    weights[n] = addWeight;
    delete[] nameProb;
    delete[] aliasStruct;
    delete[] barProb;
    delete[] blockDist;
    delete[] sortedprob;
    init(elements, weights);
    return 0;
}

int bPps::del(CountType delIdx) {
    vector<CountType> elements(n - 1);
    vector<WeightType> weights(n - 1);
    int last = 0;
    for (int i = 0; i < n; i++) {
        if (nameProb[i].first != delIdx) {
            elements[last] = nameProb[i].first;
            weights[last++] = nameProb[i].second * W;
        }
    }
    delete[] nameProb;
    delete[] aliasStruct;
    delete[] barProb;
    delete[] blockDist;
    delete[] sortedprob;
    init(elements, weights);
    return 0;
}