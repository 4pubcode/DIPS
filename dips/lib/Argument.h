#pragma once
class Argument {
   public:
    int lambda = 1;
    string setdir, pdf, resultdir, analysisdir, filename, algo,
        test = "query", graphdir, graphlabel, graphEdgepath, setpath,
        resultpath, analysispath, errpath;
    int debug = 0, times = 100, isGraph = 0, sum = 1, runtime = 0;
    CountType n = 100, updNum = 1000;
    Argument(int argc, char *argv[]);
    ~Argument();
};

Argument::Argument(int argc, char *argv[]) {
    debug = 0;
    setdir = "./data/";
    resultdir = "./result/";
    analysisdir = "./analysis/";
    graphdir = "data";
    // algo = "dss";
    // pdf = "weibull";
    // n = 100000;
    // updNum = 100000;
    // times = 10000;
    // test = "query";
    // pdf = "exp";
    // algo = "wiss";
    string param, value;
    for (int ind = 1; ind < argc; ind++) {
        if (argv[ind][0] != '-') break;

        stringstream sstr(argv[ind]);
        getline(sstr, param, '=');
        getline(sstr, value, '=');

        if (!param.compare("-pdf"))
            pdf = value;
        else if (!param.compare("-lambda"))
            lambda = stoi(value);
        else if (!param.compare("-algo"))
            algo = value;
        else if (!param.compare("-sum"))
            sum = stod(value);
        else if (!param.compare("-setdir"))
            setdir = value;
        else if (!param.compare("-debug"))
            debug = stoi(value);
        else if (!param.compare("-resultdir"))
            resultdir = value;
        else if (!param.compare("-analysisdir"))
            analysisdir = value;
        else if (!param.compare("-test"))
            test = value;
        else if (!param.compare("-times"))
            times = stoi(value);
        else if (!param.compare("-updNum"))
            updNum = stoi(value);
        else if (!param.compare("-n"))  // the number of elements
            n = stoi(value);
        // n = stoll(value.c_str());
        else if (!param.compare("-graphdir"))
            graphdir = value;
        else if (!param.compare("-graph"))
            isGraph = stoi(value);
        else if (!param.compare("-runtime"))
            runtime = stoi(value);
        else if (!param.compare("-graphlabel"))
            graphlabel = value;
        else if (!param.compare("-dist"))
            pdf = value;
    }

    if (!isGraph)
        if (pdf == "exp") {
            setpath = setdir + pdf + "_" + to_string(lambda) + "_" +
                      to_string(sum) + "_" + to_string(n) + "_" +
                      to_string(updNum);
            resultpath = resultdir + algo + "_" + pdf + "_" +
                         to_string(lambda) + "_" + to_string(sum) + "_" +
                         to_string(n) + "_" + to_string(updNum);
            analysispath = analysisdir + algo + "_" + pdf + "_" +
                           to_string(lambda) + "_" + to_string(sum) + "_" +
                           to_string(n) + "_" + to_string(updNum);
            errpath = analysisdir + "Err_" + algo + "_" + pdf + "_" +
                      to_string(lambda) + "_" + to_string(sum) + "_" +
                      to_string(n) + "_" + to_string(updNum);
        } else {
            setpath = setdir + pdf + "_" + to_string(sum) + "_" + to_string(n) +
                      "_" + to_string(updNum);
            resultpath = resultdir + algo + "_" + pdf + "_" + to_string(sum) +
                         "_" + to_string(n) + "_" + to_string(updNum);
            analysispath = analysisdir + algo + "_" + pdf + "_" +
                           to_string(sum) + "_" + to_string(n) + "_" +
                           to_string(updNum);
            errpath = analysisdir + "Err_" + algo + "_" + pdf + "_" +
                      to_string(sum) + "_" + to_string(n) + "_" +
                      to_string(updNum);
        }
    else {
        string graphAttr = graphdir + graphlabel + ".attribute";
        ifstream graphAttrIn(graphAttr.c_str());
        CountType vertexNum, edgeNum;
        string tmp;
        graphAttrIn >> tmp >> vertexNum >> tmp >> edgeNum;
        cout << "vertexNum=" << vertexNum << " edgeNum=" << edgeNum << endl;
        if (vertexNum == 0) {
            cout << "read attribute error" << endl;
            exit(-1);
        }
        n = edgeNum;
        graphEdgepath = graphdir + graphlabel + "_" + pdf + ".outWEdges_new";
        resultpath = resultdir + algo + "_" + graphlabel + "_" + pdf + "_" +
                     to_string(updNum);
        analysispath = analysisdir + "World_" + algo + "_" + graphlabel + "_" +
                       pdf + "_" + to_string(updNum);
    }
    // setpath = setdir + "normal_-1_10000000000_10000";
}

Argument::~Argument() {}
