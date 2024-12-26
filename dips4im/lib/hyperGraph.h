#pragma once
template <class T>
class HyperGraph {
   private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;
    /// _numRRsets: number of RR sets.
    size_t _numRRsets = 0;
    vector<bool> _vecVisitBool;
    Nodelist _vecVisitNode;

    size_t _hit = 0;
    double _numSamplesEval = 0;
    double _hyperedgeAvgEval = 0.0;
    /// Initialization
    void InitHypergraph() {
        _numV = (uint32_t)_graph.n;
        _numE = _graph.m;
        _FRsets = FRsets(_numV);
        _vecVisitBool = vector<bool>(_numV);
        _vecVisitNode = Nodelist(_numV);
    }

   public:
    /// _graph: reverse graph
    T &_graph;
    /// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can
    /// reach
    FRsets _FRsets;
    /// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach
    /// node i
    RRsets _RRsets;

    explicit HyperGraph(T &graph) : _graph(graph) { InitHypergraph(); }

    /// Returns the number of nodes in the graph.
    uint32_t get_nodes() { return _numV; }

    /// Returns the number of edges in the graph.
    size_t get_edges() { return _numE; }

    /// Returns the number of RR sets in the graph.
    size_t get_RR_sets_size() { return _numRRsets; }

    void RefreshHypergraph() {
        if (_RRsets.size() != 0) {
            for (auto i = _numRRsets; i--;) {
                RRset().swap(_RRsets[i]);
            }

            RRsets().swap(_RRsets);

            for (auto i = _numV; i--;) {
                FRset().swap(_FRsets[i]);
            }
        }

        _numRRsets = 0;
        _hit = 0;
    }

    /// Generate a set of n RR sets
    void BuildRRsets(size_t numSamples) {
        void (*func)(uint32_t uStart, size_t hyperIdx);

        if (numSamples > SIZE_MAX) {
            cout << "Error:R too large" << endl;
            exit(1);
        }

        auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        for (auto i = prevSize; i < numSamples; i++) {
            CountType numVisitNode =
                _graph.BuildOneRRsetSkewed(_vecVisitBool, _vecVisitNode);
            for (auto j = 0; j < numVisitNode; j++)
                _FRsets[_vecVisitNode[j]].push_back(i);
            _RRsets.push_back(RRset(_vecVisitNode.begin(),
                                    _vecVisitNode.begin() + numVisitNode));
        }

        return;
    }

    // Evaluate the influence spread of a seed set on current generated RR sets
    double CalculateInf(Nodelist &vecSeed) {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        for (auto seed : vecSeed) {
            for (auto node : _FRsets[seed]) {
                vecBoolVst[node] = true;
            }
        }

        size_t count = std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
        return 1.0 * count * _numV / _numRRsets;
    }

    // Release memory
    void ReleaseMemory() {
        RefreshHypergraph();
        vector<bool>().swap(_vecVisitBool);
        Nodelist().swap(_vecVisitNode);
        FRsets().swap(_FRsets);
    }
};
