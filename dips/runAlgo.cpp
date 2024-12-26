#include "lib/DIPS.h"
#include "lib/bPps.h"
#include "lib/dPps.h"
#include "lib/hPps.h"
#include "lib/nPps.h"
#include "lib/ppsMethod.h"
#include "lib/sPps.h"
int Logger::mode;

template <class T>
void runQuery(T &algo, vector<WeightType> &weights, CountType times,
              double &time, vector<CountType> &result) {
    auto n = weights.size();
    vector<CountType> elements(n);
    for (int i = 0; i < n; i++) {
        elements[i] = i;
    }
    algo.init(elements, weights);
    // cout << "init finish" << endl;
    Timer queryTimer;
    queryTimer.get_operation_time();
    for (int i = 0; i < times; i++) {
        algo.query(result);
    }
    time = queryTimer.get_operation_time() / double(times);
}

template <class T>
void runUpdate(T &algo, CountType n, vector<WeightType> &weights,
               CountType toAdd, vector<WeightType> &add, CountType toDel,
               vector<CountType> &del, double &time,
               vector<CountType> &result) {
    vector<CountType> elements(n);
    for (int i = 0; i < n; i++) {
        elements[i] = i;
    }
    algo.init(elements, weights);
    cout << "init finish" << endl;
    Timer updateTimer;
    updateTimer.get_operation_time();
    for (CountType j = 0; j < toAdd; j++) {
        algo.del(del[j]);
        algo.add(add[j], j + n);
    }
    time = updateTimer.get_operation_time();
    algo.query(result);
    time = time / double(toAdd) / 2;
}

int main(int argc, char *argv[]) {
    init_random_seed();
    Argument Arg(argc, argv);
    Logger::mode = Arg.debug;
    vector<WeightType> weights(Arg.n);
    CountType toAdd, toDel;
    vector<WeightType> add;
    vector<CountType> del;

    /* Load data */
    ifstream in(Arg.setpath + ".weights", ios::in | ios::binary);
    if (!in) {
        Logger::LogError("infile not found", Arg.setpath + ".weights");
        exit(-1);
    }
    in.read((char *)&weights[0], Arg.n * sizeof(weights[0]));
    // add
    if (Arg.test == "update" || Arg.test == "all") {
        toAdd = Arg.updNum;
        add.resize(toAdd);
        in.read((char *)&add[0], toAdd * sizeof(add[0]));
    }
    in.close();
    // del
    if (Arg.test == "update" || Arg.test == "all") {
        toDel = toAdd;
        del.resize(toDel);
        ifstream inDel(Arg.setpath + ".del", ios::in | ios::binary);
        if (!inDel) {
            Logger::LogError("inDelfile not found", Arg.setpath + ".del");
            exit(-1);
        }
        inDel.read((char *)&del[0], toDel * sizeof(del[0]));
        inDel.close();
    }

    /* Run query and/or update */
    string outalgoname;
    int times = Arg.times;
    vector<CountType> result;
    if (Arg.test == "query" || Arg.test == "all") {
        result.resize(0);
        double time;
        if (Arg.algo == "dips") {
            DIPS dips;
            runQuery(dips, weights, times, time, result);
            outalgoname = "dips";
        } else if (Arg.algo == "dss") {
            dPps dpps;
            runQuery(dpps, weights, times, time, result);
            outalgoname = "dpps";
        } else if (Arg.algo == "bss") {
            bPps bss;
            runQuery(bss, weights, times, time, result);
            outalgoname = "bringmann";
        } else if (Arg.algo == "nss") {
            naivePps nss;
            runQuery(nss, weights, times, time, result);
            outalgoname = "naive";
        } else if (Arg.algo == "sss") {
            staticPps sss;
            runQuery(sss, weights, times, time, result);
            outalgoname = "static";
        } else if (Arg.algo == "hss") {
            hPps hss;
            runQuery(hss, weights, times, time, result);
            outalgoname = "hybrid";
        }
        ofstream out(Arg.resultpath + ".out");
        if (!out) {
            Logger::LogError("query outfile cannot be created",
                             Arg.resultpath + ".out");
            exit(-1);
        }
        for (size_t i = 0; i < result.size(); i++) out << result[i] << endl;
        out.close();
        ofstream outtime(Arg.analysispath + ".runtime", ios::app);
        outtime << time << "," << endl;
        outtime.close();
    }
    if (Arg.test == "update" || Arg.test == "all") {
        toAdd = 1;
        result.resize(0);
        double time;
        Timer updateTimer;
        if (Arg.algo == "dips") {
            DIPS dips;
            runUpdate(dips, Arg.n, weights, toAdd, add, toDel, del, time,
                      result);
        }
        // else if (Arg.algo == "dpps")
        // {
        //     mdSubsetSampling dpps;
        //     runUpdate(dpps, Arg.n, weights, toAdd, add, toDel, del, time,
        //     result); outalgoname = "dpps";
        // }
        else if (Arg.algo == "dss") {
            dPps dpps;
            runUpdate(dpps, Arg.n, weights, toAdd, add, toDel, del, time,
                      result);
        } else if (Arg.algo == "bss") {
            bPps bss;
            runUpdate(bss, Arg.n, weights, toAdd, add, toDel, del, time,
                      result);
        } else if (Arg.algo == "nss") {
            naivePps nss;
            runUpdate(nss, Arg.n, weights, toAdd, add, toDel, del, time,
                      result);
        } else if (Arg.algo == "sss") {
            staticPps sss;
            runUpdate(sss, Arg.n, weights, toAdd, add, toDel, del, time,
                      result);
        } else if (Arg.algo == "hss") {
            hPps hss;
            runUpdate(hss, Arg.n, weights, toAdd, add, toDel, del, time,
                      result);
        }
        ofstream out(Arg.resultpath + ".out");
        if (!out) {
            Logger::LogError(
                "update outfile cannot be created",
                Arg.resultpath + outalgoname + "_" + Arg.filename + ".out");
            exit(-1);
        }
        for (size_t i = 0; i < result.size(); i++) out << result[i] << endl;
        out.close();
        ofstream outtime(Arg.analysispath + ".updtime", ios::app);
        outtime << time << "," << endl;
        outtime.close();
    }
    return 0;
}
