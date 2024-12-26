#include "lib/alg.cpp"
#include "lib/bPps.h"
#include "lib/common_graph.h"
#include "lib/hPps.h"
#include "lib/nPps.h"
#include "lib/symbols.h"
#include "lib/DIPS.h"
#include "lib/dPps.h"
// #include "lib/dss_graph.h"

template <class T>
void run(T graph, double &time, int seedSize, double eps) {
    Alg<T> tAlg(graph);
    double delta = 1.0 / graph.n;
    time = tAlg.subsimOnly(seedSize, eps, delta);
}

lookup_table *g_tables = new lookup_table(2, BASE);
vector<vector<int>> *g_bits2set = new vector<vector<int>>(256);

int main(int argc, char *argv[]) {
    Argument Arg(argc, argv);
    // Arg.algo = "bss";
    // Arg._seedsize = 100;
    // Arg.pdf = "exp";
    // Arg.filelabel = "orkut_exp";
    init_random_seed();
    int seedSize = Arg._seedsize;
    double time;
    string s;
    CountType upd_num = 1000;

    std::cout << "seedSize k=" << seedSize << std::endl;
    std::cout << "starting..." << std::endl;
    string fname = Arg.filelabel + '_' + Arg.algo + '_' +
                   to_string(seedSize) + ".updtime";
    ofstream oupdtime(Arg.analysispath + fname, ios::app);
    if (Arg.algo == "dips") {
        for (int i = 0; i < 256; i++) {
            auto tmp = i;
            for (int j = 0; j < 8; j++) {
                if (tmp % 2 == 1) {
                    (*g_bits2set)[i].push_back(j);
                }
                tmp >>= 1;
            }
        }
        common_originGraph<DIPS> g(Arg.filedir, Arg.filelabel,

                                               Arg.pdf);
        // s = "wiss_";
        run(g, time, seedSize, Arg._eps);
        oupdtime << g.update(upd_num) << endl;
    } else if (Arg.algo == "mdss") {
        common_originGraph<dPps> g(Arg.filedir, Arg.filelabel,

                                                    Arg.pdf);
        // dss_originGraph g(Arg.filedir, Arg.filelabel, Arg.pdf);
        // s = "mdss_";
        run(g, time, seedSize, Arg._eps);
        oupdtime << g.update(upd_num) << endl;
    } else if (Arg.algo == "nss") {
        common_originGraph<nPps> g(Arg.filedir, Arg.filelabel,

                                                  Arg.pdf);
        // nss_originGraph g(Arg.filedir, Arg.filelabel, Arg.pdf);
        // s = "nss_";
        run(g, time, seedSize, Arg._eps);
        oupdtime << g.update(upd_num) << endl;
    } else if (Arg.algo == "hss") {
        common_originGraph<hPps> g(Arg.filedir, Arg.filelabel,

                                                   Arg.pdf);
        // hss_originGraph g(Arg.filedir, Arg.filelabel, Arg.pdf);
        // s = "hss_";
        run(g, time, seedSize, Arg._eps);
        oupdtime << g.update(upd_num) << endl;
    } else if (Arg.algo == "bss") {
        common_originGraph<bPps> g(Arg.filedir,
                                                      Arg.filelabel,

                                                      Arg.pdf);
        // bss_originGraph g(Arg.filedir, Arg.filelabel, Arg.pdf);
        // s = "bss_";
        run(g, time, seedSize, Arg._eps);
        oupdtime << g.update(upd_num) << endl;
    } else if (Arg.algo == "dss") {
        common_originGraph<dPps> g(Arg.filedir, Arg.filelabel,

                                                    Arg.pdf);
        // dss_originGraph g(Arg.filedir, Arg.filelabel, Arg.pdf);
        // s = "dss_";
        run(g, time, seedSize, Arg._eps);
        oupdtime << g.update(upd_num) << endl;
    }
    oupdtime.close();
    /* Output running time to analysis dir */
    // s += Arg.filelabel + '_' + Arg.pdf + '_' + to_string(seedSize) + ".imtime";
    s = Arg.filelabel + '_' + Arg.algo + '_' +
                   to_string(seedSize) + ".imtime";
    ofstream outqtime(Arg.analysispath + s, ios::app);
    outqtime << time << "," << endl;
    outqtime.close();
    cout << "write to " << s << endl;

    return 0;
}
