#include "ppsMethod.h"
class nPps : public ppsMethod {
   private:
    CountType n;
    element *nameProb;
    // CountType *dict;
    unordered_map<CountType, CountType> dict;
    CountType newSize;
    WeightType W;

   public:
    nPps(/* args */);
    ~nPps();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &results);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
};

nPps::nPps(/* args */) {}

nPps::~nPps() {}

int nPps::query(vector<CountType> &results) {
    results.resize(0);
    if (n == 0) return 0;
    for (size_t i = 0; i < n; i++) {
        // if (dsfmt_gv_genrand_open_open()*W < nameProb[i].second) {
        if (dsfmt_gv_genrand_open_close() * W <= nameProb[i].second) {
            results.push_back(nameProb[i].first);
        }
    }
    return 0;
}

int nPps::init(vector<CountType> elements, vector<WeightType> weights) {
    n = weights.size();
    if (n == 0) return 0;
    if (n > CountTypeMAX / 2)
        newSize = CountTypeMAX;
    else
        newSize = n * 2;
    // dict = new CountType[newSize];
    nameProb = new element[newSize];
    W = 0;
    for (CountType i = 0; i < n; i++) {
        W += weights[i];
    }
    for (CountType i = 0; i < n; i++) {
        nameProb[i].first = elements[i];
        nameProb[i].second = weights[i];
        dict[elements[i]] = i;
    }
}


int nPps::add(WeightType addWeight, CountType idx) {
    // vector<CountType> elements(n+1);
    // vector<WeightType> weights(n+1);
    // auto newW = W+addWeight;
    // auto ratio =  W/newW;
    // for (int i=0; i<n; i++) {
    //     elements[i] = nameProb[i].first;
    //     weights[i] = nameProb[i].second*ratio;
    // }
    // elements[n] = idx;
    // weights[n] = addWeight/newW;
    // delete [] nameProb;
    // // delete [] dict;
    // init(elements, weights);
    n++;
    if (n > newSize) {
        cout << "error, overflow" << endl;
        exit(-1);
    }
    nameProb[n - 1].second = addWeight;
    nameProb[n - 1].first = idx;
    W += addWeight;
    dict[idx] = n - 1;
}

/*delIdx: the name of the element */
int nPps::del(CountType delIdx) {
    // vector<CountType> elements(n-1);
    // vector<WeightType> weights(n-1);
    // auto newW = W-nameProb[dict[delIdx]].second;
    // auto ratio =  W/newW;
    // int last = 0;
    // for (int i=0; i<n; i++) {
    //     if (i != dict[delIdx]) {
    //         elements[last] = nameProb[i].first;
    //         weights[last++] = nameProb[i].second*ratio;
    //     }
    // }
    // delete [] nameProb;
    // // delete [] dict;
    // init(elements, weights);
    n--;
    CountType idx = dict[delIdx];  // the index of delIdx in the namaProb
    W -= nameProb[idx].second;
    nameProb[idx].first = nameProb[n].first;
    nameProb[idx].second = nameProb[n].second;
    dict[nameProb[n].first] = idx;  // the n-th element is moved to idx
    dict[delIdx] = -1;
}