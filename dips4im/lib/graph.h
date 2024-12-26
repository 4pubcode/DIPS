#pragma once 
#include <chrono>
class originGraph
{
public:
    uint n, m;
    string filedir, filelabel;
    vector<node> *neighborList;
    unordered_map<uint, uint> *adjList;
    string dist;
    uint *outSizeList;
    originGraph() {}

    int generateDist(vector<WeightType> &result, string dist, CountType nums, WeightType sum)
    {
        result.resize(nums);
        default_random_engine generator(time(NULL));
        if (dist == "weibull")
        {
            double min_value = 1e-8;
            double sum = 0.0;
            for (size_t j = 0; j < nums; j++)
            {
                // random number from (0, 10)
                double a = dsfmt_gv_genrand_open_open() * 10;
                double b = dsfmt_gv_genrand_open_open() * 10;
                weibull_distribution<double> distribution(a, b);
                auto weight = distribution(generator);
                result[j] = weight;
                sum += weight;
            }
            if (sum < 0)
            {
                cout << "******************error" << endl;
            }
            for (size_t j = 0; j < nums; j++)
            {
                auto weight = result[j] / sum;
                result[j] = (weight > min_value) ? weight : min_value;
            }
            // sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
        }
        else if (dist == "exp")
        {
            double min_value = 1e-8;
            double sum = 0.0;
            for (size_t j = 0; j < nums; j++)
            {
                // lambda = 1
                auto weight = -log(1.0 - dsfmt_gv_genrand_open_open());
                result[j] = weight;
                sum += weight;
            }
            for (size_t j = 0; j < nums; j++)
            {
                double weight = result[j] / sum;
                result[j] = (weight > min_value) ? weight : min_value;
            }
            // sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
        }
        else if (dist == "normal")
        {
            double min_value = 1e-8;
            normal_distribution<double> distribution(0, 1);
            double sum = 0.0;
            for (size_t j = 0; j < nums; j++)
            {
                // lambda = 1
                auto weight = distribution(generator);
                result[j] = weight;
            }
            double max = -100000000, min = 100000000;
            for (size_t j = 0; j < nums; j++)
            {
                if (result[j] > max)
                    max = result[j];
                if (result[j] < min)
                    min = result[j];
            }
            for (size_t j = 0; j < nums; j++)
            {
                result[j] = result[j] - min;
                if (result[j] < min_value)
                    result[j] = min_value;
                if (result[j] < 0)
                {
                    cout << "error min" << endl;
                    exit(-1);
                }
                sum += result[j];
            }
            for (size_t j = 0; j < nums; j++)
            {
                double weight = result[j] / sum;
                result[j] = (weight > min_value) ? weight : min_value;
            }
            // sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
        }
    }

    originGraph(const string &_filedir, const string &_filelabel, const string &_dist)
    {
        filedir = _filedir + _filelabel + "_weighted/";
        filelabel = _filelabel;
        dist = _dist;
        string neiNode, neiWeight, neiNum, graphAttr;
        neiNode = filedir + filelabel + ".outEdges";
        neiWeight = filedir + filelabel + "_" + dist + ".outWEdges_new";
        neiNum = filedir + filelabel + ".outPtr";
        graphAttr = filedir + filelabel + ".attribute";
        cout << "FilePath: " << graphAttr.c_str() << endl;
        cout << "Read graph attributes..." << endl;
        string tmp;
        ifstream graphAttrIn(graphAttr.c_str());
        if (!graphAttrIn)
        {
            uint m = 0;
            string filename = filedir + filelabel + ".txt";
            ifstream infile;
            infile.open(filename);
            if (!infile)
            {
                cout << "ERROR: unable to open source dataset: " << filename << endl;
                exit(-1);
            }

            cout << "Read the original txt file..." << endl;
            uint from;
            uint to;
            double weight;
            uint n = 0;
            while (infile >> from >> to >> weight)
            {
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
            for (uint i = 0; i < n; i++)
            {
                indegree[i] = 0;
                outdegree[i] = 0;
            }
            // read graph and get degree info

            while (infile >> from >> to >> weight)
            {
                outdegree[from]++;
                indegree[to]++;
            }

            uint **inAdjList = new uint *[n];
            uint **outAdjList = new uint *[n];

            double **inAdjWList = new double *[n];
            double **outAdjWList = new double *[n];

            uint *pointer_in = new uint[n];
            uint *pointer_out = new uint[n];
            for (uint i = 0; i < n; i++)
            {
                inAdjList[i] = new uint[indegree[i]];
                outAdjList[i] = new uint[outdegree[i]];

                inAdjWList[i] = new double[indegree[i]];
                outAdjWList[i] = new double[outdegree[i]];

                pointer_in[i] = 0;
                pointer_out[i] = 0;
            }
            infile.clear();
            infile.seekg(0);

            while (infile >> from >> to >> weight)
            {
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
            for (uint i = 0; i < n; i++)
            {
                outid += outdegree[i];
                coutPL[i + 1] = outid;
                for (uint j = 0; j < outdegree[i]; j++)
                {
                    coutEL[out_curnum] = outAdjList[i][j];
                    coutWEL[out_curnum] = outAdjWList[i][j];
                    out_curnum += 1;
                }
            }
            ofstream cout_attr(filedir + filelabel + ".attribute");
            cout_attr << "n " << n << "\n";
            cout_attr << "m " << m << "\n";
            cout_attr.close();
            ofstream foutEL(filedir + filelabel + ".outEdges", ios::out | ios::binary);
            ofstream foutPL(filedir + filelabel + ".outPtr", ios::out | ios::binary);
            ofstream foutWEL(filedir + filelabel + "_" + dist + ".outWEdges_new", ios::out | ios::binary);
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
        if (n == 0)
        {
            cout << "read attribute error" << endl;
            exit(-1);
        }
        graphAttrIn.close();
        cout << "Read graph ..." << endl;
        ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
        ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
        ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
        if (!neiWeightIn)
        {
            uint outSize;
            uint outSizeSum = 0, preOutSizeSum = 0;
            neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
            cout << "first outSizesum:" << outSizeSum << endl;
            ofstream neiWeightOut(neiWeight.c_str(), ios::out | ios::binary);
            for (uint i = 0; i < n; i++)
            {
                neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
                outSize = outSizeSum - preOutSizeSum;
                preOutSizeSum = outSizeSum;
                vector<WeightType> result;
                generateDist(result, dist, outSize, 1);
                neiWeightOut.write((char *)&result[0], sizeof(result[0]) * outSize);
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
        neighborList = new vector<node>[n];
        adjList = new unordered_map<uint, uint>[n];
        outSizeList = new uint[n];
        for (uint i = 0; i < n; i++)
        {
            neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
            outSize = outSizeSum - preOutSizeSum;
            preOutSizeSum = outSizeSum;
            for (uint j = 0; j < outSize; j++)
            {
                uint id;
                double w;
                neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
                neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
                neighborList[i].push_back(node(id, w));
                adjList[i][id] = j;
            }
            outSizeList[i] = outSize;
        }
        neiNumIn.close();
        neiWeightIn.close();
        neiNodeIn.close();
        cout << "finish graph init" << endl;
    }

    double update(CountType add_num)
    {
        uint *sarr = new uint[add_num];
        uint *tarr = new uint[add_num];
        double *warr = new double[add_num];
        int opnum = 0;
        ifstream opfile(filedir + filelabel + "_" + dist + '_' + to_string(add_num) + ".op_new", ios::in);
        if (!opfile)
        {
            cout << "Generate add edges.." << endl;
            opfile.open(filedir + "/" + filelabel + '_' + to_string(add_num) + ".op", ios::in);
            if (!opfile)
            {
                getAddEdge(add_num);
                opfile.open(filedir + "/" + filelabel + '_' + to_string(add_num) + ".op", ios::in);
            }
            while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
            {
                opnum++;
            }
            for (int i = 0; i < add_num; i++)
            {
                if (adjList[sarr[i]].find(tarr[i]) == adjList[sarr[i]].end())
                {
                    cout << "delete a nonexistent neighbor" << endl;
                    exit(-2);
                }
                warr[i] = neighborList[sarr[i]][adjList[sarr[i]][tarr[i]]].w;
            }
            opfile.close();
            ofstream out;
            out.open(filedir + filelabel + "_" + dist + '_' + to_string(add_num) + ".op_new", ios::out);
            for (int i = 0; i < add_num; i++)
            {
                out << sarr[i] << " " << tarr[i] << " " << warr[i] << endl;
            }
            out.close();
            opfile.open(filedir + filelabel + "_" + dist + '_' + to_string(add_num) + ".op_new", ios::in);
        }
        opnum = 0;
        while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
        {
            opnum++;
        }

        cout << "the first addition" << sarr[0] << " " << tarr[0] << endl;
        opfile.close();
        // delete
        auto begin = chrono::high_resolution_clock::now();
        for (int i = 0; i < add_num; i++)
        {
            uint s = sarr[i];
            uint t = tarr[i];
            uint neiidx = adjList[s][t];
            uint outsize = outSizeList[s];
            if (neiidx == outsize - 1)
                neighborList[s].pop_back();
            else
            {
                auto tmp = neighborList[s].end() - 1;
                neighborList[s][neiidx].id = tmp->id;
                neighborList[s][neiidx].w = tmp->w;
                adjList[s][tmp->id] = neiidx;
                neighborList[s].pop_back();
            }
            adjList[s].erase(t);
            outSizeList[s]--;
        }
        for (int i = 0; i < add_num; i++)
        {
            uint s = sarr[i];
            neighborList[s].push_back(node{tarr[i], warr[i]});
            outSizeList[s]++;
            adjList[s][tarr[i]] = neighborList[s].size() - 1;
        }
        return (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / (add_num * 2);
    }

    ~originGraph()
    {
    }

    void getAddEdge(long num)
    {
        cout << "generate edges to add" << endl;
        pair<int, double> *aliasD = new pair<int, double>[n];
        unordered_map<int, vector<int>> hasAdded;
        for (int i = 0; i < n; i++)
        {
            aliasD[i] = make_pair(i, outSizeList[i]);
        }
        Alias alias = Alias(aliasD, n);
        ofstream output(filedir + filelabel + ".op", ios::out);
        int i = num;
        while (i)
        {
            int nodeidx = alias.generateRandom(dsfmt_gv_genrand_close_open(), dsfmt_gv_genrand_close_open());
            int outSize = outSizeList[nodeidx];
            while (outSize <= 1)
            {
                nodeidx = alias.generateRandom(dsfmt_gv_genrand_close_open(), dsfmt_gv_genrand_close_open());
                outSize = outSizeList[nodeidx];
            }
            int neiidx = floor(dsfmt_gv_genrand_close_open() * outSize);
            auto tmp = neighborList[nodeidx][neiidx];
            if (hasAdded[nodeidx].size() > 0)
            {
                int flag = 0;
                for (auto j : hasAdded[nodeidx])
                    if (j == tmp.id)
                    {
                        flag = 1;
                        break;
                    }
                if (flag == 1)
                    continue;
            }
            output << nodeidx << " " << tmp.id << " " << tmp.w << endl;
            hasAdded[nodeidx].push_back(tmp.id);
            i--;
        }
        delete[] aliasD;
        output.close();
        cout << "finish edges added" << endl;
    }
};
