#include <fstream>
#include <iostream>
#include <random>
#include "SFMT/dSFMT/dSFMT.c"
#include "SFMT/dSFMT/dSFMT.h"
using namespace std;

int main(void) {
    string ds = "orkut";
    mt19937 generator(time(NULL));
    ifstream infile;
    string filename = "data/" + ds + ".txt";
    infile.open(filename, std::ios::in);
    if (!infile) {
        cout << "ERROR: unable to open source dataset: " << filename << endl;
        exit(-1);
    }
    int from;
    int to;
    int n = 5;
    string ofname = "data/" + ds + "_exp.txt";
    ofstream outfile;
    // double min_w = 100000;
    // vector<double> weights;
    // int rows = 117185087;
    // for (int i=0; i<rows; i++) {
    //     double a = dsfmt_gv_genrand_open_open() * 10;
    //     double b = dsfmt_gv_genrand_open_open() * 10;
    //     weibull_distribution<double> distribution(a, b);
    //     auto weight = distribution(generator);
    //     weights.push_back(weight);
    //     if (weight < min_w) {
    //         min_w = weight;
    //     }
    // }

    int ind = 0;
    outfile.open(ofname, std::ios::out);
    exponential_distribution<double> distribution(1);
    if (!outfile) {
        cout << "ERROR: unable to open"<< endl;
        exit(-1);
    }
    while (infile >> from >> to) {
        auto weight = 1 + distribution(generator);
        // auto weight = 1 + rand();
        outfile << from << " " << to << " " << weight << endl;
        // outfile << from << " " << to << " " << 1+weights[ind++]-min_w << endl;
    }
}