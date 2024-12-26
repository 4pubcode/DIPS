#include "lib/symbols.h"
int Logger::mode;

void init_random_seed(uint32_t seed) {
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(seed));
}

int main(int argc, char *argv[]) {
    Argument Arg(argc, argv);
    init_random_seed(Arg.runtime);
    Logger::mode = Arg.debug;
    CountType totGen = Arg.n + Arg.updNum;  // Arg.n /100 for udpate
    double *weights = new double[totGen];
    double min_value = (1e-8 < 1.0 / Arg.n) ? 1e-8 : 1.0 / Arg.n;
    mt19937 generator(time(NULL));

    if (Arg.pdf == "weibull") {
        double a = 5;
        double b = 1;
        weibull_distribution<double> distribution(a, b);
        for (size_t i = 0; i < totGen; i++) {
            auto weight = distribution(generator);
            weights[i] = weight;
        }
    } else if (Arg.pdf == "exp") {
        for (size_t i = 0; i < totGen; i++) {
            exponential_distribution<> distribution(0.5);
            auto weight = distribution(generator);
            weights[i] = weight;
        }
    } else if (Arg.pdf == "gamma") {
        for (size_t i = 0; i < totGen; i++) {
            gamma_distribution<double> distribution(2.0, 2.0);
            auto weight = distribution(generator);
            weights[i] = weight;
        }
    } else if (Arg.pdf == "normal") {
        for (size_t i = 0; i < totGen; i++) {
            normal_distribution<double> distribution(0, 10);
            auto weight = distribution(generator);
            weights[i] = weight;
        }
    } else if (Arg.pdf == "halfnormal") {
        for (size_t i = 0; i < totGen; i++) {
            normal_distribution<double> distribution(0, 10);
            auto weight = fabs(distribution(generator));
            weights[i] = weight;
        }
    } else if (Arg.pdf == "lognormal") {
        for (size_t i = 0; i < totGen; i++) {
            lognormal_distribution<double> distribution(0, sqrt(log(2)));
            auto weight = distribution(generator);
            weights[i] = weight;
        }
    } else if (Arg.pdf == "uniform") {
        for (size_t i = 0; i < totGen; i++) {
            weights[i] = dsfmt_gv_genrand_open_open();
        }
    }

    double max = -100000000, min = 100000000;
    for (size_t i = 0; i < Arg.n; i++) {
        if (weights[i] > max) max = weights[i];
        if (weights[i] < min) min = weights[i];
    }
    double gap = max - min;
    for (size_t i = 0; i < totGen; i++) {
        weights[i] = weights[i] - min + 1;
        if (weights[i] < 0) {
            Logger::LogError("error min");
            exit(-1);
        }
    }

    ofstream out(Arg.setpath + ".weights", ios::out | ios::binary);
    if (!out) {
        Logger::LogError("open out weights file error:", Arg.filename);
        exit(-1);
    }
    out.write((char *)&weights[0], sizeof(weights[0]) * totGen);
    out.close();

    Logger::LogDebug("The first 100 numbers:\n");
    for (size_t i = 0; i < 100; i++) Logger::LogDebug(i, weights[i]);

    // generate random number for delete
    // if (Arg.n >= 1000)
    // {
    CountType toDelNum = Arg.updNum;
    vector<CountType> toDel;
    vector<int> got(Arg.n + Arg.updNum, 0);
    for (CountType i = 0; i < toDelNum; i++) {
        CountType tmp = floor(dsfmt_gv_genrand_open_open() * (Arg.n + i));
        while (got[tmp] != 0)
            tmp = floor(dsfmt_gv_genrand_open_open() * (Arg.n + i));
        toDel.push_back(tmp);
        got[tmp] = 1;
    }
    out.open(Arg.setpath + ".del", ios::out | ios::binary);
    if (!out) {
        Logger::LogError("open out del file error", Arg.filename);
        exit(-1);
    }
    out.write((char *)&toDel[0], sizeof(toDel[0]) * toDelNum);
    out.close();
    // }

    return 0;
}