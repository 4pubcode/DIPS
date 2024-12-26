#pragma once
#include <string>
using namespace std;
class Argument {
   public:
    string filedir, filelabel, resultpath, analysispath, algo, pdf;
    int debug = 0, times = 100, updNum = 1000, _seedsize = 5;
    double _eps = 0.1;
    Argument(int argc, char *argv[]);
    ~Argument();
};

Argument::Argument(int argc, char *argv[]) {
    debug = 0;
    resultpath = "./result/";
    analysispath = "./analysis/";
    filedir = "./data/";
    filelabel = "orkut";
    algo = "dss";
    pdf = "weibull";
    string param, value;
    for (int ind = 1; ind < argc; ind++) {
        if (argv[ind][0] != '-') break;
        stringstream sstr(argv[ind]);
        getline(sstr, param, '=');
        getline(sstr, value, '=');
        if (!param.compare("-filedir"))
            filedir = value;
        else if (!param.compare("-filelabel"))
            filelabel = value;
        else if (!param.compare("-algo"))
            algo = value;
        else if (!param.compare("-debug"))
            debug = stoi(value);
        else if (!param.compare("-resultpath"))
            resultpath = value;
        else if (!param.compare("-pdf"))
            pdf = value;
        else if (!param.compare("-analysispath"))
            analysispath = value;
        else if (!param.compare("-times"))
            times = stoi(value);
        else if (!param.compare("-updNum"))
            updNum = stoi(value);
        else if (!param.compare("-seedsize"))
            _seedsize = stoi(value);
    }
}

Argument::~Argument() {}
