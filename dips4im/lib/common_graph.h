#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "ppsMethod.h"
#include "symbols.h"

using namespace std;
template <class SSM>
class common_originGraph {
   public:
    uint n, m;
    string filedir, filelabel;
    node **neighborList;
    // vector<vector<node>> neighborList;
    // unordered_map<uint, uint> *adjList;
    shared_ptr<unordered_map<uint, uint>[]> adjList;
    // SSM *ssm;
    shared_ptr<SSM[]> ssm;
    // vector<SSM> ssm;
    string dist;
    shared_ptr<uint[]> ngbSize;
    // uint *ngbSize;
    // CountType *newSize;
    common_originGraph() {}
    CountType BuildOneRRsetSkewed(vector<bool> &_vecVisitBool,
                                  vector<uint> &_vecVisitNode) {
        CountType uStart = floor(dsfmt_gv_genrand_open_open() * n);
        CountType numVisitNode = 0, currIdx = 0;
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        while (currIdx < numVisitNode) {
            const auto curNode = _vecVisitNode[currIdx++];
            vector<CountType> sample_results;
            auto numNeighbor = ngbSize[curNode];
            if (numNeighbor > 0) {
                ssm[curNode].query(sample_results);
                for (auto id : sample_results) {
                    if (!_vecVisitBool[id]) {
                        _vecVisitNode[numVisitNode++] = id;
                        _vecVisitBool[id] = true;
                    }
                }
            }
            // for (auto res:sample_results) {
            //     auto newNode = neighborList[curNode][res];
            //      _vecVisitNode[numVisitNode++] = newNode.id;
            //      _vecVisitBool[newNode.id] = true;
            // }
            // TODO: doing sample
            //  auto numNeighbor = outSizeList[tempNode];
            //  for (CountType j = 0; j < numNeighbor; j++)
            //  {
            //      double prob = dsfmt_gv_genrand_close_open();
            //      auto newNode = neighborList[tempNode][j];
            //      if (prob > newNode.w || _vecVisitBool[newNode.id])
            //      {
            //          continue;
            //      }
            //      _vecVisitNode[numVisitNode++] = newNode.id;
            //      _vecVisitBool[newNode.id] = true;
            //  }
        }

        for (uint i = 0; i < numVisitNode; i++)
            _vecVisitBool[_vecVisitNode[i]] = false;
        return numVisitNode;
    }

    WeightType generateDist(vector<WeightType> &result, const string &dist,
                            CountType nums) {
        result.resize(nums);
        WeightType sum = 0;
        default_random_engine generator(time(NULL));
        if (dist == "weibull") {
            for (size_t j = 0; j < nums; j++) {
                // random number from (0, 10)
                double a = dsfmt_gv_genrand_open_open() * 10;
                double b = dsfmt_gv_genrand_open_open() * 10;
                weibull_distribution<double> distribution(a, b);
                auto weight = 1 + distribution(generator);
                result[j] = weight;
                sum += weight;
            }
        } else if (dist == "exp") {
            for (size_t j = 0; j < nums; j++) {
                // lambda = 1
                auto weight = 1 - log(1.0 - dsfmt_gv_genrand_open_open());
                result[j] = weight;
                sum += weight;
            }
        } else if (dist == "normal") {
            normal_distribution<double> distribution(0, 1);
            for (size_t j = 0; j < nums; j++) {
                // lambda = 1
                auto weight = 1 + distribution(generator);
                result[j] = weight;
                sum += weight;
            }
        }
    }

    common_originGraph(const string &_filedir, const string &_filelabel,
                       const string &_dist) {
        // filedir = _filedir + _filelabel + "_weighted/";
        filedir = _filedir;
        filelabel = _filelabel;
        dist = _dist;
        string neiNode, neiWeight, neiNum, graphAttr;
        neiNode = filedir + filelabel + ".outEdges";
        neiWeight = filedir + filelabel + "_" + dist + ".outWEdges";
        neiNum = filedir + filelabel + ".outPtr";
        graphAttr = filedir + filelabel + ".attribute";
        cout << "FilePath: " << graphAttr.c_str() << endl;
        cout << "Read graph attributes..." << endl;
        string tmp;

        ifstream graphAttrIn(graphAttr.c_str());
        if (!graphAttrIn) {  // first time processing a dataset
            uint m = 0;
            string filename = filedir + filelabel + ".txt";
            ifstream infile;
            infile.open(filename, std::ios::in);
            if (!infile) {
                cout << "ERROR: unable to open source dataset: " << filename
                     << endl;
                exit(-1);
            }

            cout << "Read the original txt file..." << endl;
            uint from;
            uint to;
            double weight;
            uint n = 0;
            while (infile >> from >> to >> weight) {
                if (from > n)
                    n = from;
                else if (to > n)
                    n = to;
            }
            n++;
            infile.clear();
            infile.seekg(0);
            uint *indegree = new uint[n];
            uint *outdegree = new uint[n];
            for (uint i = 0; i < n; i++) {
                indegree[i] = 0;
                outdegree[i] = 0;
            }
            // read graph and get degree info

            while (infile >> from >> to >> weight) {
                outdegree[from]++;
                indegree[to]++;
            }

            uint **inAdjList = new uint *[n];
            uint **outAdjList = new uint *[n];

            double **inAdjWList = new double *[n];
            double **outAdjWList = new double *[n];

            uint *pointer_in = new uint[n];
            uint *pointer_out = new uint[n];
            for (uint i = 0; i < n; i++) {
                inAdjList[i] = new uint[indegree[i]];
                outAdjList[i] = new uint[outdegree[i]];

                inAdjWList[i] = new double[indegree[i]];
                outAdjWList[i] = new double[outdegree[i]];

                pointer_in[i] = 0;
                pointer_out[i] = 0;
            }

            // second pass
            infile.clear();
            infile.seekg(0);

            while (infile >> from >> to >> weight) {
                outAdjList[from][pointer_out[from]] = to;
                outAdjWList[from][pointer_out[from]] = weight;
                pointer_out[from]++;
                inAdjList[to][pointer_in[to]] = from;
                inAdjWList[to][pointer_in[to]] = weight;
                pointer_in[to]++;
                m++;
            }
            infile.close();
            delete pointer_in;
            delete pointer_out;

            cout << "Write to csr file..." << endl;
            uint *coutEL = new uint[m];
            uint *coutPL = new uint[n + 1];
            double *coutWEL = new double[m];
            coutPL[0] = 0;
            uint outid = 0;
            uint out_curnum = 0;
            for (uint i = 0; i < n; i++) {
                outid += outdegree[i];
                coutPL[i + 1] = outid;
                for (uint j = 0; j < outdegree[i]; j++) {
                    coutEL[out_curnum] = outAdjList[i][j];
                    coutWEL[out_curnum] = outAdjWList[i][j];
                    out_curnum += 1;
                }
            }
            ofstream cout_attr(filedir + filelabel + ".attribute");
            cout_attr << "n " << n << "\n";
            cout_attr << "m " << m << "\n";
            cout_attr.close();
            ofstream foutEL(filedir + filelabel + ".outEdges",
                            ios::out | ios::binary);
            ofstream foutPL(filedir + filelabel + ".outPtr",
                            ios::out | ios::binary);
            // ofstream foutWEL(
            //     filedir + filelabel + "_" + dist + ".outWEdges_new",
            //     ios::out | ios::binary);
            ofstream foutWEL(filedir + filelabel + "_" + dist + ".outWEdges",
                             ios::out | ios::binary);
            foutEL.write((char *)&coutEL[0], sizeof(coutEL[0]) * m);
            foutPL.write((char *)&coutPL[0], sizeof(coutPL[0]) * (n + 1));
            foutWEL.write((char *)&coutWEL[0], sizeof(coutWEL[0]) * m);

            foutEL.close();
            foutPL.close();
            foutWEL.close();

            delete[] coutPL;
            delete[] coutEL;
            delete[] coutWEL;
            graphAttrIn.close();
            graphAttrIn.open(graphAttr.c_str());
        }
        graphAttrIn >> tmp >> n >> tmp >> m;
        cout << "n=" << n << endl;
        if (n == 0) {
            cout << "read attribute error" << endl;
            exit(-1);
        }
        graphAttrIn.close();

        cout << "Read graph ..." << endl;
        ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
        ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
        ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
        if (!neiWeightIn) {
            uint outSize;
            uint outSizeSum = 0, preOutSizeSum = 0;
            neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
            cout << "first outSizesum:" << outSizeSum << endl;
            ofstream neiWeightOut(neiWeight.c_str(), ios::out | ios::binary);
            for (uint i = 0; i < n; i++) {
                neiNumIn.read(reinterpret_cast<char *>(&outSizeSum),
                              sizeof(uint));
                outSize = outSizeSum - preOutSizeSum;
                preOutSizeSum = outSizeSum;
                vector<WeightType> result;
                auto sum = generateDist(result, dist, outSize);
                neiWeightOut.write((char *)&result[0],
                                   sizeof(result[0]) * outSize);
            }
            neiWeightOut.close();
            neiWeightIn.open(neiWeight.c_str(), ios::in | ios::binary);
            neiNumIn.close();
            neiNumIn.open(neiNum.c_str(), ios::in | ios::binary);
        }

        uint outSize;
        uint outSizeSum = 0, preOutSizeSum = 0;
        neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
        cout << "first outSizesum:" << outSizeSum << endl;

        neighborList = new node *[n];
        // neighborList = vector<vector<node>>(n);
        // newSize = new CountType[n];
        // adjList = new unordered_map<uint, uint>[n];
        adjList = shared_ptr<unordered_map<uint, uint>[]>(
            new unordered_map<uint, uint>[n]);
        // ssm = new SSM[n];
        ssm = shared_ptr<SSM[]>(new SSM[n]);
        // ssm = vector<SSM>(n);
        // ngbSize = new uint[n];
        ngbSize = shared_ptr<uint[]>(new uint[n]);
        for (uint i = 0; i < n; i++) {
            neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
            outSize = outSizeSum - preOutSizeSum;
            preOutSizeSum = outSizeSum;
            CountType newSizetmp;
            if (outSize > CountTypeMAX / 2)
                newSizetmp = CountTypeMAX;
            else
                newSizetmp = outSize * 2;
            neighborList[i] = new node[newSizetmp];
            // neighborList[i] = vector<node>(newSizetmp);
            vector<WeightType> weights;
            vector<CountType> elements;
            for (uint j = 0; j < outSize; j++) {
                uint id;
                double w;
                neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
                neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
                neighborList[i][j] = node(id, w);
                elements.push_back(id);
                weights.push_back(w);
                adjList[i][id] = j;
            }

            ssm[i].init(elements, weights);
            ngbSize[i] = outSize;
        }
        neiNumIn.close();
        neiWeightIn.close();
        neiNodeIn.close();
        cout << "finish graph init" << endl;
    }

    double update(CountType add_num) {
        uint *sarr = new uint[add_num];
        uint *tarr = new uint[add_num];
        double *warr = new double[add_num];
        int opnum = 0;
        ifstream opfile(filedir + filelabel + "_" + dist + '_' +
                            to_string(add_num) + ".op_new",
                        ios::in);
        if (!opfile) {
            cout << "Generate add edges.." << endl;
            opfile.open(filedir + filelabel + '_' + to_string(add_num) + ".op",
                        ios::in);
            if (!opfile) {
                getAddEdge(add_num);
                opfile.open(
                    filedir + filelabel + '_' + to_string(add_num) + ".op",
                    ios::in);
            }
            while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum]) {
                opnum++;
            }
            for (int i = 0; i < add_num; i++) {
                auto s = sarr[i];
                auto t = tarr[i];
                if (adjList[s].find(t) == adjList[s].end()) {
                    cout << "delete a nonexistent neighbor" << endl;
                    exit(-2);
                }
                warr[i] = neighborList[s][adjList[s][t]].w;
            }
            opfile.close();

            ofstream out;
            out.open(filedir + filelabel + "_" + dist + '_' +
                         to_string(add_num) + ".op_new",
                     ios::out);
            for (int i = 0; i < add_num; i++) {
                out << sarr[i] << " " << tarr[i] << " " << warr[i] << endl;
            }
            out.close();

            opfile.open(filedir + filelabel + "_" + dist + '_' +
                            to_string(add_num) + ".op_new",
                        ios::in);
        }
        opnum = 0;
        while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum]) {
            opnum++;
        }
        opfile.close();
        cout << "the first addition: " << sarr[0] << " " << tarr[0] << endl;

        // delete
        auto begin = chrono::high_resolution_clock::now();
        for (int i = 0; i < add_num; i++) {
            uint s = sarr[i];
            uint t = tarr[i];
            uint neiidx = adjList[s][t];
            uint outsize = ngbSize[s];
            if (neiidx == outsize - 1)
                neighborList[s][neiidx].id = -1;
            else {
                auto tmp = outsize - 1;
                neighborList[s][neiidx].id = neighborList[s][tmp].id;
                neighborList[s][neiidx].w = neighborList[s][tmp].w;
                adjList[s][neighborList[s][tmp].id] = neiidx;
                neighborList[s][tmp].id = -1;
            }
            adjList[s].erase(t);
            ngbSize[s]--;
            ssm[s].del(t);
        }
        for (int i = 0; i < add_num; i++) {
            uint s = sarr[i];
            uint t = tarr[i];
            uint outsize = ngbSize[s];
            // if (outsize + 1 > newSize[s]) {
            //     cout << "error, overflow" << endl;
            //     exit(-1);
            // }
            neighborList[s][outsize] = node{tarr[i], warr[i]};
            ngbSize[s]++;
            adjList[s][tarr[i]] = outsize;
            ssm[s].add(warr[i], t);
        }
        auto time = (chrono::duration_cast<chrono::nanoseconds>(
                         chrono::high_resolution_clock::now() - begin)
                         .count() /
                     1000000000.0) /
                    (add_num * 2);
        delete[] sarr;
        delete[] tarr;
        delete[] warr;
        return time;
    }

    ~common_originGraph() {
        // delete []ssm;
        // delete[] adjList;
        // delete[] newSize;
    }

    void getAddEdge(long num) {
        cout << "generate edges to add" << endl;
        pair<int, double> *aliasD = new pair<int, double>[n];
        unordered_map<int, vector<int>> hasAdded;
        for (int i = 0; i < n; i++) {
            aliasD[i] = make_pair(i, ngbSize[i]);
        }
        Alias alias = Alias(aliasD, n);
        ofstream ofile(filedir + filelabel + '_' + to_string(num) + ".op",
                       ios::out);

        while (num) {
            int nodeidx;
            int outSize = 0;
            while (outSize <= 1) {
                nodeidx = alias.generateRandom(dsfmt_gv_genrand_close_open(),
                                               dsfmt_gv_genrand_close_open());
                outSize = ngbSize[nodeidx];
            }

            int neiidx = floor(dsfmt_gv_genrand_close_open() * outSize);
            auto ngb = neighborList[nodeidx][neiidx];
            int flag = 0;
            if (hasAdded[nodeidx].size() > 0) {
                for (auto v : hasAdded[nodeidx]) {
                    if (v == ngb.id) {
                        flag = 1;
                        break;
                    }
                }
            }
            if (flag == 1) continue;
            ofile << nodeidx << " " << ngb.id << " " << ngb.w << endl;
            hasAdded[nodeidx].push_back(ngb.id);
            num--;
        }
        delete[] aliasD;
        ofile.close();
        cout << "finish edges added" << endl;
    }
};

// #pragma once
// #include <fstream>
// #include <iostream>

// #include "subsetSamplingMethod.h"
// #include "symbols.h"

// #define BASE 4
// using namespace std;
// template <class SSM>
// class common_originGraph {
//    public:
//     uint n, m;
//     string filedir, filelabel;
//     node **neighborList;
//     // vector<vector<node>> neighborList;
//     unordered_map<uint, uint> *adjList;
//     SSM *ssm;
//     // vector<SSM> ssm;
//     string dist;
//     uint *ngbSize;
//     // CountType *newSize;
//     common_originGraph() {}
//     CountType BuildOneRRsetSkewed(vector<bool> &_vecVisitBool,
//                                   vector<uint> &_vecVisitNode) {
//         CountType uStart = floor(dsfmt_gv_genrand_open_open() * n);
//         CountType numVisitNode = 0, currIdx = 0;
//         _vecVisitNode[numVisitNode++] = uStart;
//         _vecVisitBool[uStart] = true;
//         while (currIdx < numVisitNode) {
//             const auto curNode = _vecVisitNode[currIdx++];
//             vector<CountType> sample_results;
//             auto numNeighbor = ngbSize[curNode];
//             if (numNeighbor > 0) {
//                 ssm[curNode].query(sample_results);
//                 for (auto id : sample_results) {
//                     if (!_vecVisitBool[id]) {
//                         _vecVisitNode[numVisitNode++] = id;
//                         _vecVisitBool[id] = true;
//                     }
//                 }
//             }
//             // for (auto res:sample_results) {
//             //     auto newNode = neighborList[curNode][res];
//             //      _vecVisitNode[numVisitNode++] = newNode.id;
//             //      _vecVisitBool[newNode.id] = true;
//             // }
//             // TODO: doing sample
//             //  auto numNeighbor = outSizeList[tempNode];
//             //  for (CountType j = 0; j < numNeighbor; j++)
//             //  {
//             //      double prob = dsfmt_gv_genrand_close_open();
//             //      auto newNode = neighborList[tempNode][j];
//             //      if (prob > newNode.w || _vecVisitBool[newNode.id])
//             //      {
//             //          continue;
//             //      }
//             //      _vecVisitNode[numVisitNode++] = newNode.id;
//             //      _vecVisitBool[newNode.id] = true;
//             //  }
//         }

//         for (uint i = 0; i < numVisitNode; i++)
//             _vecVisitBool[_vecVisitNode[i]] = false;
//         return numVisitNode;
//     }

//     WeightType generateDist(vector<WeightType> &result, const string &dist,
//                             CountType nums) {
//         result.resize(nums);
//         WeightType sum = 0;
//         default_random_engine generator(time(NULL));
//         if (dist == "weibull") {
//             for (size_t j = 0; j < nums; j++) {
//                 // random number from (0, 10)
//                 double a = dsfmt_gv_genrand_open_open() * 10;
//                 double b = dsfmt_gv_genrand_open_open() * 10;
//                 weibull_distribution<double> distribution(a, b);
//                 auto weight = 1 + distribution(generator);
//                 result[j] = weight;
//                 sum += weight;
//             }
//         } else if (dist == "exp") {
//             for (size_t j = 0; j < nums; j++) {
//                 // lambda = 1
//                 auto weight = 1 - log(1.0 - dsfmt_gv_genrand_open_open());
//                 result[j] = weight;
//                 sum += weight;
//             }
//         } else if (dist == "normal") {
//             normal_distribution<double> distribution(0, 1);
//             for (size_t j = 0; j < nums; j++) {
//                 // lambda = 1
//                 auto weight = 1 + distribution(generator);
//                 result[j] = weight;
//                 sum += weight;
//             }
//         }
//     }

//     common_originGraph(const string &_filedir, const string &_filelabel,
//                        const string &_dist) {
//         // filedir = _filedir + _filelabel + "_weighted/";
//         filedir = _filedir;
//         filelabel = _filelabel;
//         dist = _dist;
//         string neiNode, neiWeight, neiNum, graphAttr;
//         neiNode = filedir + filelabel + ".outEdges";
//         neiWeight = filedir + filelabel + "_" + dist + ".outWEdges";
//         neiNum = filedir + filelabel + ".outPtr";
//         graphAttr = filedir + filelabel + ".attribute";
//         cout << "FilePath: " << graphAttr.c_str() << endl;
//         cout << "Read graph attributes..." << endl;
//         string tmp;

//         ifstream graphAttrIn(graphAttr.c_str());
//         if (!graphAttrIn) {  // first time processing a dataset
//             uint m = 0;
//             string filename = filedir + filelabel + ".txt";
//             ifstream infile;
//             infile.open(filename, std::ios::in);
//             if (!infile) {
//                 cout << "ERROR: unable to open source dataset: " << filename
//                      << endl;
//                 exit(-1);
//             }

//             cout << "Read the original txt file..." << endl;
//             uint from;
//             uint to;
//             double weight;
//             uint n = 0;
//             while (infile >> from >> to >> weight) {
//                 if (from > n)
//                     n = from;
//                 else if (to > n)
//                     n = to;
//             }
//             n++;
//             infile.clear();
//             infile.seekg(0);
//             uint *indegree = new uint[n];
//             uint *outdegree = new uint[n];
//             for (uint i = 0; i < n; i++) {
//                 indegree[i] = 0;
//                 outdegree[i] = 0;
//             }
//             // read graph and get degree info

//             while (infile >> from >> to >> weight) {
//                 outdegree[from]++;
//                 indegree[to]++;
//             }

//             uint **inAdjList = new uint *[n];
//             uint **outAdjList = new uint *[n];

//             double **inAdjWList = new double *[n];
//             double **outAdjWList = new double *[n];

//             uint *pointer_in = new uint[n];
//             uint *pointer_out = new uint[n];
//             for (uint i = 0; i < n; i++) {
//                 inAdjList[i] = new uint[indegree[i]];
//                 outAdjList[i] = new uint[outdegree[i]];

//                 inAdjWList[i] = new double[indegree[i]];
//                 outAdjWList[i] = new double[outdegree[i]];

//                 pointer_in[i] = 0;
//                 pointer_out[i] = 0;
//             }

//             // second pass
//             infile.clear();
//             infile.seekg(0);

//             while (infile >> from >> to >> weight) {
//                 outAdjList[from][pointer_out[from]] = to;
//                 outAdjWList[from][pointer_out[from]] = weight;
//                 pointer_out[from]++;
//                 inAdjList[to][pointer_in[to]] = from;
//                 inAdjWList[to][pointer_in[to]] = weight;
//                 pointer_in[to]++;
//                 m++;
//             }
//             infile.close();
//             delete pointer_in;
//             delete pointer_out;

//             cout << "Write to csr file..." << endl;
//             uint *coutEL = new uint[m];
//             uint *coutPL = new uint[n + 1];
//             double *coutWEL = new double[m];
//             coutPL[0] = 0;
//             uint outid = 0;
//             uint out_curnum = 0;
//             for (uint i = 0; i < n; i++) {
//                 outid += outdegree[i];
//                 coutPL[i + 1] = outid;
//                 for (uint j = 0; j < outdegree[i]; j++) {
//                     coutEL[out_curnum] = outAdjList[i][j];
//                     coutWEL[out_curnum] = outAdjWList[i][j];
//                     out_curnum += 1;
//                 }
//             }
//             ofstream cout_attr(filedir + filelabel + ".attribute");
//             cout_attr << "n " << n << "\n";
//             cout_attr << "m " << m << "\n";
//             cout_attr.close();
//             ofstream foutEL(filedir + filelabel + ".outEdges",
//                             ios::out | ios::binary);
//             ofstream foutPL(filedir + filelabel + ".outPtr",
//                             ios::out | ios::binary);
//             // ofstream foutWEL(
//             //     filedir + filelabel + "_" + dist + ".outWEdges_new",
//             //     ios::out | ios::binary);
//             ofstream foutWEL(
//                 filedir + filelabel + "_" + dist + ".outWEdges",
//                 ios::out | ios::binary);
//             foutEL.write((char *)&coutEL[0], sizeof(coutEL[0]) * m);
//             foutPL.write((char *)&coutPL[0], sizeof(coutPL[0]) * (n + 1));
//             foutWEL.write((char *)&coutWEL[0], sizeof(coutWEL[0]) * m);

//             foutEL.close();
//             foutPL.close();
//             foutWEL.close();

//             delete[] coutPL;
//             delete[] coutEL;
//             delete[] coutWEL;
//             graphAttrIn.close();
//             graphAttrIn.open(graphAttr.c_str());
//         }
//         graphAttrIn >> tmp >> n >> tmp >> m;
//         cout << "n=" << n << endl;
//         if (n == 0) {
//             cout << "read attribute error" << endl;
//             exit(-1);
//         }
//         graphAttrIn.close();

//         cout << "Read graph ..." << endl;
//         ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
//         ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
//         ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
//         if (!neiWeightIn) {
//             uint outSize;
//             uint outSizeSum = 0, preOutSizeSum = 0;
//             neiNumIn.read(reinterpret_cast<char *>(&outSizeSum),
//             sizeof(uint)); cout << "first outSizesum:" << outSizeSum << endl;
//             ofstream neiWeightOut(neiWeight.c_str(), ios::out | ios::binary);
//             for (uint i = 0; i < n; i++) {
//                 neiNumIn.read(reinterpret_cast<char *>(&outSizeSum),
//                               sizeof(uint));
//                 outSize = outSizeSum - preOutSizeSum;
//                 preOutSizeSum = outSizeSum;
//                 vector<WeightType> result;
//                 auto sum = generateDist(result, dist, outSize);
//                 neiWeightOut.write((char *)&result[0],
//                                    sizeof(result[0]) * outSize);
//             }
//             neiWeightOut.close();
//             neiWeightIn.open(neiWeight.c_str(), ios::in | ios::binary);
//             neiNumIn.close();
//             neiNumIn.open(neiNum.c_str(), ios::in | ios::binary);
//         }

//         uint outSize;
//         uint outSizeSum = 0, preOutSizeSum = 0;
//         neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
//         cout << "first outSizesum:" << outSizeSum << endl;

//         neighborList = new node *[n];
//         // neighborList = vector<vector<node>>(n);
//         // newSize = new CountType[n];
//         adjList = new unordered_map<uint, uint>[n];
//         ssm = new SSM[n];
//         // ssm = vector<SSM>(n);
//         ngbSize = new uint[n];
//         for (uint i = 0; i < n; i++) {
//             neiNumIn.read(reinterpret_cast<char *>(&outSizeSum),
//             sizeof(uint)); outSize = outSizeSum - preOutSizeSum;
//             preOutSizeSum = outSizeSum;
//             CountType newSizetmp;
//             if (outSize > CountTypeMAX / 2)
//                 newSizetmp = CountTypeMAX;
//             else
//                 newSizetmp = outSize * 2;
//             neighborList[i] = new node[newSizetmp];
//             // neighborList[i] = vector<node>(newSizetmp);
//             vector<WeightType> weights;
//             vector<CountType> elements;
//             for (uint j = 0; j < outSize; j++) {
//                 uint id;
//                 double w;
//                 neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
//                 neiWeightIn.read(reinterpret_cast<char *>(&w),
//                 sizeof(double)); neighborList[i][j] = node(id, w);
//                 elements.push_back(id);
//                 weights.push_back(w);
//                 adjList[i][id] = j;
//             }

//             ssm[i].init(elements, weights);
//             ngbSize[i] = outSize;
//         }
//         neiNumIn.close();
//         neiWeightIn.close();
//         neiNodeIn.close();
//         cout << "finish graph init" << endl;
//     }

//     double update(CountType add_num) {
//         uint *sarr = new uint[add_num];
//         uint *tarr = new uint[add_num];
//         double *warr = new double[add_num];
//         int opnum = 0;
//         ifstream opfile(filedir + filelabel + "_" + dist + '_' +
//                             to_string(add_num) + ".op_new",
//                         ios::in);
//         if (!opfile) {
//             cout << "Generate add edges.." << endl;
//             opfile.open(
//                 filedir + filelabel + '_' + to_string(add_num) + ".op",
//                 ios::in);
//             if (!opfile) {
//                 getAddEdge(add_num);
//                 opfile.open(filedir + filelabel + '_' +
//                                 to_string(add_num) + ".op",
//                             ios::in);
//             }
//             while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum]) {
//                 opnum++;
//             }
//             for (int i = 0; i < add_num; i++) {
//                 auto s = sarr[i];
//                 auto t = tarr[i];
//                 if (adjList[s].find(t) == adjList[s].end()) {
//                     cout << "delete a nonexistent neighbor" << endl;
//                     exit(-2);
//                 }
//                 warr[i] = neighborList[s][adjList[s][t]].w;
//             }
//             opfile.close();

//             ofstream out;
//             out.open(filedir + filelabel + "_" + dist + '_' +
//                          to_string(add_num) + ".op_new",
//                      ios::out);
//             for (int i = 0; i < add_num; i++) {
//                 out << sarr[i] << " " << tarr[i] << " " << warr[i] << endl;
//             }
//             out.close();

//             opfile.open(filedir + filelabel + "_" + dist + '_' +
//                             to_string(add_num) + ".op_new",
//                         ios::in);
//         }
//         opnum = 0;
//         while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum]) {
//             opnum++;
//         }
//         opfile.close();
//         cout << "the first addition: " << sarr[0] << " " << tarr[0] << endl;

//         // delete
//         auto begin = chrono::high_resolution_clock::now();
//         for (int i = 0; i < add_num; i++) {
//             uint s = sarr[i];
//             uint t = tarr[i];
//             uint neiidx = adjList[s][t];
//             uint outsize = ngbSize[s];
//             if (neiidx == outsize - 1)
//                 neighborList[s][neiidx].id = -1;
//             else {
//                 auto tmp = outsize - 1;
//                 neighborList[s][neiidx].id = neighborList[s][tmp].id;
//                 neighborList[s][neiidx].w = neighborList[s][tmp].w;
//                 adjList[s][neighborList[s][tmp].id] = neiidx;
//                 neighborList[s][tmp].id = -1;
//             }
//             adjList[s].erase(t);
//             ngbSize[s]--;
//             ssm[s].del(t);
//         }
//         for (int i = 0; i < add_num; i++) {
//             uint s = sarr[i];
//             uint t = tarr[i];
//             uint outsize = ngbSize[s];
//             // if (outsize + 1 > newSize[s]) {
//             //     cout << "error, overflow" << endl;
//             //     exit(-1);
//             // }
//             neighborList[s][outsize] = node{tarr[i], warr[i]};
//             ngbSize[s]++;
//             adjList[s][tarr[i]] = outsize;
//             ssm[s].add(warr[i],t);
//         }
//         return (chrono::duration_cast<chrono::nanoseconds>(
//                     chrono::high_resolution_clock::now() - begin)
//                     .count() /
//                 1000000000.0) /
//                (add_num * 2);

//         delete []sarr;
//         delete []tarr;
//         delete []warr;
//     }

//     ~common_originGraph() {
//         delete []ssm;
//         // delete[] adjList;
//         // delete[] newSize;
//     }

//     void getAddEdge(long num) {
//         cout << "generate edges to add" << endl;
//         pair<int, double> *aliasD = new pair<int, double>[n];
//         unordered_map<int, vector<int>> hasAdded;
//         for (int i = 0; i < n; i++) {
//             aliasD[i] = make_pair(i, ngbSize[i]);
//         }
//         Alias alias = Alias(aliasD, n);
//         ofstream ofile(filedir + filelabel + '_' + to_string(num) + ".op",
//         ios::out);

//         while (num) {
//             int nodeidx;
//             int outSize = 0;
//             while (outSize <= 1) {
//                 nodeidx = alias.generateRandom(dsfmt_gv_genrand_close_open(),
//                                                dsfmt_gv_genrand_close_open());
//                 outSize = ngbSize[nodeidx];
//             }

//             int neiidx = floor(dsfmt_gv_genrand_close_open() * outSize);
//             auto ngb = neighborList[nodeidx][neiidx];
//             int flag = 0;
//             if (hasAdded[nodeidx].size() > 0) {
//                 for (auto v : hasAdded[nodeidx]) {
//                     if (v == ngb.id) {
//                         flag = 1;
//                         break;
//                     }
//                 }
//             }
//             if (flag == 1) continue;
//             ofile << nodeidx << " " << ngb.id << " " << ngb.w << endl;
//             hasAdded[nodeidx].push_back(ngb.id);
//             num--;
//         }
//         delete[] aliasD;
//         ofile.close();
//         cout << "finish edges added" << endl;
//     }
// };
