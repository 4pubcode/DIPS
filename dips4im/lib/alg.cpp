#include "alg.h"
// code from OPIM
template <class T>
double Alg<T>::MaxCoverVanilla(const int targetSize) {
    // optimization with minimum upper bound among all rounds [Default].
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;) {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    // degMap: map degree to the nodes with this degree
    RRsets degMap(maxDeg + 1);

    for (auto i = _numV; i--;) {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--)  // Enusre deg > 0
    {
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;) {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];  //

            if (deg > currDeg) {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true) {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--) {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound) {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--) {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--) {
                    const auto currDegBound =
                        coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound) {
                        flag = false;
                    } else if (vecBound[0] < currDegBound) {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound) {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--) {
                        const auto currDegBound =
                            coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound) {
                            flag = false;
                        } else if (vecBound[0] < currDegBound) {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(),
                                               size_t(0)) +
                                    sumInf) *
                             _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize) {
                // Top-k influential nodes constructed
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // cout << ">>>[greedy-lazy] influence: " << finalInf << ",
                // min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx]) {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx]) {
                    if (coverage[nodeIdx] == 0)
                        continue;  // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV;  // All RR sets are covered.
}

template <class T>
double Alg<T>::MaxCoverTopK(const int targetSize) {
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;) {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg +
                  1);  // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;) {
        // if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }

    Nodelist sortedNode(_numV);    // sortedNode: record the sorted nodes in
                                   // ascending order of degree
    Nodelist nodePosition(_numV);  // nodePosition: record the position of each
                                   // node in the sortedNode
    Nodelist degreePosition(
        maxDeg +
        2);  // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;

    for (auto &nodes : degMap) {
        degreePosition[idxDegree + 1] =
            degreePosition[idxDegree] + (uint32_t)nodes.size();
        idxDegree++;

        for (auto &node : nodes) {
            nodePosition[node] = idxSort;
            sortedNode[idxSort++] = node;
        }
    }

    // check if an edge is removed
    vector<bool> edgeMark(_numRRsets, false);
    // record the total of top-k marginal gains
    size_t sumTopk = 0;

    for (auto deg = maxDeg + 1; deg--;) {
        if (degreePosition[deg] <= _numV - targetSize) {
            sumTopk += deg * (degreePosition[deg + 1] - (_numV - targetSize));
            break;
        }

        sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }

    _boundMin = 1.0 * sumTopk;
    _vecSeed.clear();
    size_t sumInf = 0;

    /*
     * sortedNode: position -> node
     * nodePosition: node -> position
     * degreePosition: degree -> position (start position of this degree)
     * coverage: node -> degree
     * e.g., swap the position of a node with the start position of its degree
     * swap(sortedNode[nodePosition[node]],
     * sortedNode[degreePosition[coverage[node]]])
     */
    for (auto k = targetSize; k--;) {
        const auto seed = sortedNode.back();
        sortedNode.pop_back();
        const auto newNumV = sortedNode.size();
        sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
        sumInf += coverage[seed];
        _vecSeed.push_back(seed);
        coverage[seed] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[seed]) {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx]) {
                if (coverage[nodeIdx] == 0)
                    continue;  // This node is seed, skip

                const auto currPos =
                    nodePosition[nodeIdx];               // The current position
                const auto currDeg = coverage[nodeIdx];  // The current degree
                const auto startPos =
                    degreePosition[currDeg];  // The start position of this
                                              // degree
                const auto startNode =
                    sortedNode[startPos];  // The node with the start position
                // Swap this node to the start position with the same degree,
                // and update their positions in nodePosition
                swap(sortedNode[currPos], sortedNode[startPos]);
                nodePosition[nodeIdx] = startPos;
                nodePosition[startNode] = currPos;
                // Increase the start position of this degree by 1, and decrease
                // the degree of this node by 1
                degreePosition[currDeg]++;
                coverage[nodeIdx]--;

                // If the start position of this degree is in top-k, reduce topk
                // by 1
                if (startPos >= newNumV - targetSize) sumTopk--;
            }
        }

        _boundLast = 1.0 * (sumInf + sumTopk);

        if (_boundMin > _boundLast) _boundMin = _boundLast;
    }

    _boundMin *= 1.0 * _numV / _numRRsets;
    _boundLast *= 1.0 * _numV / _numRRsets;
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    cout << "  >>>[greedy-topk] influence: " << finalInf
         << ", min-bound: " << _boundMin << ", last-bound: " << _boundLast
         << '\n';
    return finalInf;
}

template <class T>
double Alg<T>::MaxCover(const int targetSize) {
    if (targetSize >= 1000) return MaxCoverTopK(targetSize);

    return MaxCoverVanilla(targetSize);
}

template <class T>
double Alg<T>::subsimOnly(const int targetSize, const double epsilon,
                          const double delta) {
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta =
        sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) /
                                targetSize / pow2(epsilon)) +
                         1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    cout << "a1: " << a1 << ", a2: " << a2 << endl;

    cout << endl;
    for (auto idx = 1; idx <= numIter; idx++) {
        const auto numR = numRbase << (idx - 1);
        cout << "Iteration: " << idx << " RR set: " << numR << endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR);      // R1
        _hyperGraphVldt.BuildRRsets(numR);  // R2
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCover(targetSize);  // 从R1选S*得到的影响力
        time2 += timerSubsim.get_operation_time();
        const auto infVldt =
            _hyperGraphVldt.CalculateInf(_vecSeed);  // S*在R2的影响力

        const auto degVldt =
            infVldt * _numRRsets / _numV;  // S*在R2覆盖的RR set数量
        auto upperBound =
            _boundMin;  // R1估计的最大影响力的upper bound，应该和infVldt

        const auto upperDegOPT =
            upperBound * _numRRsets /
            _numV;  // 估计最大影响力对应覆盖的R1中RR set数量
        const auto lowerSelect =
            pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) -
            a1 / 18.0;  //
        const auto upperOPT =
            pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;
        cout << "SUBSIM (" << idx << "/" << numIter
             << ") approx. (max-cover): " << currApprox << endl
             << "lower bound: " << (lowerSelect * _numV / _numRRsets)
             << ", upperBound: " << upperBound << endl
             << "S* cover R1: " << infSelf * _numRRsets / _numV << "/"
             << _numRRsets << " = " << infSelf * _numRRsets / _numV / _numRRsets
             << ", S* cover R2: " << degVldt << "/" << _numRRsets << " = "
             << degVldt / _numRRsets << endl
             << "S* Inf in R1: " << infSelf << ", S* Inf in R2:" << infVldt
             << endl
             << " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets
             << '\n';
        // double avgSize = _hyperGraph.HyperedgeAvg();

        if (currApprox >= approx - epsilon) {
            Logger::LogInfo("approximation", currApprox);
            Logger::LogInfo("running_time", timerSubsim.get_total_time());
            Logger::LogInfo("influence", infVldt);
            Logger::LogInfo("influence_original", infSelf);
            Logger::LogInfo("RR_sets_size", _numRRsets * 2);
            cout << "==>Influence via R2: " << infVldt << ", time: " << '\n';
            cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2
                 << '\n';
            return time1 + time2;
        }
    }

    return 0;
}
