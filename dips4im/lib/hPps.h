#include "math.h"
#include "ppsMethod.h"
class hPps : public ppsMethod {
   private:
    CountType n;
    pair<CountType, WeightType> *nameProb, *X, *Y;
    double pm, sum, p_mu, p_thres;
    CountType Xsize = 0, Ysize = 0;
    CountType complement;
    CountType newSize;
    // CountType *dict;
    unordered_map<CountType, CountType> dict;
    WeightType W;
    int _add(WeightType addProb, CountType idx);
    int _del(CountType delIdx);

   public:
    hPps(/* args */);
    ~hPps();
    int init(vector<CountType> elements, vector<WeightType> weights);
    int query(vector<CountType> &results);
    int add(WeightType addWeight, CountType idx);
    int del(CountType delIdx);
};

hPps::hPps(/* args */) {}

hPps::~hPps() {}

int hPps::query(vector<CountType> &results) {
    mt19937 generator(time(NULL));
    results.resize(0);
    if (n == 0) return 0;
    for (CountType i = 0; i < Ysize; i++) {
        if (dsfmt_gv_genrand_open_close() < Y[i].second) {
            if (complement == 0) results.push_back(Y[i].first);
        } else {
            if (complement == 1) results.push_back(Y[i].first);
        }
    }
    CountType nX = Xsize;
    binomial_distribution<> distribution(nX, pm);
    auto k = distribution(generator);
    for (CountType i = 0; i < k; i++) {
        double r = dsfmt_gv_genrand_open_close();
        CountType idx = floor(r * (nX - i));
        if (idx >= nX - i) idx = nX - i - 1;
        swap(X[idx + i], X[i]);
        if (dsfmt_gv_genrand_open_close() < X[i].second / pm) {
            if (complement == 0) results.push_back(X[i].first);
        } else {
            if (complement == 1) results.push_back(X[i].first);
        }
    }
    if (complement == 1)
        for (CountType i = k; i < Xsize; i++) results.push_back(X[i].first);
    return 0;
}

int hPps::init(vector<CountType> elements, vector<WeightType> weights) {
    n = weights.size();
    if (n == 0) return 0;
    if (n > CountTypeMAX / 2)
        newSize = CountTypeMAX;
    else
        newSize = n * 2;
    nameProb = new pair<CountType, WeightType>[newSize];
    X = new pair<CountType, WeightType>[newSize];
    Y = new pair<CountType, WeightType>[newSize];
    // dict = new CountType[newSize];
    Xsize = 0, Ysize = 0;
    sum = 0;
    WeightType tmpPm = 0;
    W = 0;
    for (CountType i = 0; i < n; i++) {
        W += weights[i];
    }
    for (CountType i = 0; i < n; i++) {
        weights[i] /= W;
    }
    for (CountType i = 0; i < n; i++) {
        sum += weights[i];
        nameProb[i] = make_pair(elements[i], weights[i]);
        dict[elements[i]] = i;
        if (weights[i] > tmpPm) tmpPm = weights[i];
    }
    pm = 0;
    p_mu = sum / n;
    if (1 - tmpPm > tmpPm) {
        complement = 0;
        p_thres = sqrt(p_mu);
        for (CountType i = 0; i < n; i++) {
            if (nameProb[i].second < p_thres) {
                X[Xsize++] = nameProb[i];
                if (nameProb[i].second > pm) pm = nameProb[i].second;
            } else {
                Y[Ysize++] = nameProb[i];
            }
        }
    } else {
        complement = 1;
        p_mu = 1 - p_mu;
        p_thres = sqrt(p_mu);
        for (CountType i = 0; i < n; i++) {
            if (1 - nameProb[i].second < p_thres) {
                X[Xsize++] =
                    make_pair(nameProb[i].first, 1 - nameProb[i].second);
                if (1 - nameProb[i].second > pm) pm = 1 - nameProb[i].second;
            } else {
                Y[Ysize++] =
                    make_pair(nameProb[i].first, 1 - nameProb[i].second);
            }
        }
    }
    // cout << "=========init results==========" << endl;
    // cout << "pm:" << pm << endl;
    // cout << "p_mu:" << p_mu << endl;
    // cout << "p_thres:" << p_thres << endl;
    // cout << "total size:" << n << endl;
    // cout << "XSize:" << Xsize << endl;
    // cout << "YSize:" << Ysize << endl;
    // cout << "==============================" << endl;
}

int hPps::_add(WeightType addProb, CountType idx) {
    n++;
    if (n > newSize) {
        cout << "error, overflow" << endl;
        exit(-1);
    }

    sum += addProb;
    nameProb[n - 1] = make_pair(idx, addProb);
    dict[idx] = n - 1;

    Xsize = 0, Ysize = 0;
    WeightType tmpPm = 0;
    for (CountType i = 0; i < n; i++) {
        if (nameProb[i].second > tmpPm) tmpPm = nameProb[i].second;
    }
    pm = 0;
    p_mu = sum / n;
    if (1 - tmpPm > tmpPm) {
        complement = 0;
        p_thres = sqrt(p_mu);
        for (CountType i = 0; i < n; i++) {
            if (nameProb[i].second < p_thres) {
                X[Xsize++] = nameProb[i];
                if (nameProb[i].second > pm) pm = nameProb[i].second;
            } else {
                Y[Ysize++] = nameProb[i];
            }
        }
    } else {
        complement = 1;
        p_mu = 1 - p_mu;
        p_thres = sqrt(p_mu);
        for (CountType i = 0; i < n; i++) {
            if (1 - nameProb[i].second < p_thres) {
                X[Xsize++] =
                    make_pair(nameProb[i].first, 1 - nameProb[i].second);
                if (1 - nameProb[i].second > pm) pm = 1 - nameProb[i].second;
            } else {
                Y[Ysize++] =
                    make_pair(nameProb[i].first, 1 - nameProb[i].second);
            }
        }
    }
}

int hPps::_del(CountType delIdx) {
    n--;
    sum -= nameProb[dict[delIdx]].second;
    nameProb[dict[delIdx]] = nameProb[n];
    dict[nameProb[n].first] = dict[delIdx];
    dict[delIdx] = -1;

    Xsize = 0, Ysize = 0;
    WeightType tmpPm = 0;
    for (CountType i = 0; i < n; i++) {
        if (nameProb[i].second > tmpPm) tmpPm = nameProb[i].second;
    }
    pm = 0;
    p_mu = sum / n;
    if (1 - tmpPm > tmpPm) {
        complement = 0;
        p_thres = sqrt(p_mu);
        for (CountType i = 0; i < n; i++) {
            if (nameProb[i].second < p_thres) {
                X[Xsize++] = nameProb[i];
                if (nameProb[i].second > pm) pm = nameProb[i].second;
            } else {
                Y[Ysize++] = nameProb[i];
            }
        }
    } else {
        complement = 1;
        p_mu = 1 - p_mu;
        p_thres = sqrt(p_mu);
        for (CountType i = 0; i < n; i++) {
            if (1 - nameProb[i].second < p_thres) {
                X[Xsize++] =
                    make_pair(nameProb[i].first, 1 - nameProb[i].second);
                if (1 - nameProb[i].second > pm) pm = 1 - nameProb[i].second;
            } else {
                Y[Ysize++] =
                    make_pair(nameProb[i].first, 1 - nameProb[i].second);
            }
        }
    }
    return 1;
}

int hPps::add(WeightType addWeight, CountType idx) {
    vector<CountType> elements(n+1);
    vector<WeightType> weights(n+1);
    for (int i=0; i<n; i++) {
        elements[i] = nameProb[i].first;
        weights[i] = nameProb[i].second*W;
    }
    elements[n] = idx;
    weights[n] = addWeight;
    delete [] nameProb;
    delete [] X;
    delete [] Y;
    // delete [] dict;
    init(elements, weights);
    // double oldW = W;
    // W += addWeight;
    // double ratio = oldW/W;
    // for (int i=0; i<n; i++) {
    //     auto newProb =  nameProb[i].second * ratio;
    //     auto idx = nameProb[i].first;
    //     _del(idx);
    //     _add(newProb, idx);
    // }
    // _add(addWeight/W, idx);
    return 1;
}

int hPps::del(CountType delIdx) {
    vector<CountType> elements(n-1);
    vector<WeightType> weights(n-1);
    int last = 0;
    for (int i=0; i<n; i++) {
        if (i != dict[delIdx]) {
            elements[last] = nameProb[i].first;
            weights[last++] = nameProb[i].second*W;
        }
    }
    delete [] nameProb;
    delete [] X;
    delete [] Y;
    // delete [] dict;
    init(elements, weights);
    // double oldW = W;
    // auto delWeight = nameProb[dict[delIdx]].second * oldW;
    // W -= delWeight;
    // double ratio = oldW/W;
    // for (int i=0; i<n; i++) {
    //     auto newProb =  nameProb[i].second * ratio;
    //     auto idx = nameProb[i].first;
    //     _del(idx);
    //     _add(newProb, idx);
    // }
    // _del(delIdx);
    return 1;
}

