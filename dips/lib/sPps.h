class staticPps : public ppsMethod {
   private:
    CountType n;
    // vector<element> nameProb;
    // vector<CountType> nameIdx;
    element *nameProb;
    CountType *nameIdx;
    WeightType W;

   public:
    staticPps(/* args */);
    ~staticPps();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &results);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
};

staticPps::staticPps(/* args */) {}

staticPps::~staticPps() {}

bool mycomp(element i, element j) { return i.second > j.second; }

int staticPps::init(vector<CountType> elements, vector<WeightType> weights) {
    n = weights.size();
    if (n == 0) return 0;
    CountType newSize;
    if (n > INT32_MAX / 2)
        newSize = INT32_MAX;
    else
        newSize = n * 2;
    nameProb = new element[newSize];
    nameIdx = new CountType[newSize];
    W = 0;
    for (CountType i = 0; i < n; i++) {
        W += weights[i];
    }
    for (CountType i = 0; i < n; i++) {
        weights[i] /= W;
        nameProb[i] = make_pair(elements[i], weights[i]);
    }
    sort(nameProb, nameProb + n, mycomp);
    for (CountType i = 0; i < n; i++) {
        nameIdx[nameProb[i].first] = i;
    }
    return 1;
}

int staticPps::query(vector<CountType> &results) {
    results.resize(0);
    if (n == 0) return 0;
    size_t i = 0;
    while (i < n) {
        if (nameProb[i].second > 0.1) {
            if (dsfmt_gv_genrand_open_close() < MY_C * nameProb[i].second)
                results.push_back(nameProb[i].first);
            i++;
            continue;
        }
        double log_prob = log2(1 - nameProb[i].second);
        double uni_prob = dsfmt_gv_genrand_open_close();
        size_t next = i + floor(log2(uni_prob) / log_prob);
        // geometric_distribution<> distribution(nameProb[i].second);
        // size_t next = i + distribution(generator);
        if (next >= n) break;
        if (dsfmt_gv_genrand_open_close() <
            nameProb[next].second / nameProb[i].second) {
            results.push_back(nameProb[next].first);
        }
        i = next + 1;
    }
    return 0;
}

int staticPps::add(WeightType addWeight, CountType idx) {
    vector<CountType> elements(n + 1);
    vector<WeightType> weights(n + 1);
    for (int i = 0; i < n; i++) {
        elements[i] = nameProb[i].first;
        weights[i] = nameProb[i].second * W;
    }
    elements[n] = idx;
    weights[n] = addWeight;
    delete[] nameProb;
    delete[] nameIdx;
    init(elements, weights);
    // n++;
    // CountType i;
    // for (i = n - 2; i >= 0 && i < n; i--)
    // {
    //     if (nameProb[i].second < addProb)
    //     {
    //         nameProb[i + 1] = nameProb[i];
    //         nameIdx[nameProb[i + 1].first] = i + 1;
    //     }
    //     else
    //         break;
    // }
    // nameProb[i + 1] = make_pair(idx, addProb);
    // nameIdx[idx] = i + 1;
    return 1;
}

int staticPps::del(CountType delIdx) {
    vector<CountType> elements(n - 1);
    vector<WeightType> weights(n - 1);
    int last = 0;
    for (int i = 0; i < n; i++) {
        if (i != nameIdx[delIdx]) {
            elements[last] = nameProb[i].first;
            weights[last++] = nameProb[i].second * W;
        }
    }
    delete[] nameProb;
    delete[] nameIdx;
    init(elements, weights);
    // CountType idx = nameIdx[delIdx];
    // for (CountType i = idx; i < n; i++)
    // {
    //     nameProb[i] = nameProb[i + 1];
    //     nameIdx[nameProb[i].first] = i;
    // }
    // n--;
    return 1;
}