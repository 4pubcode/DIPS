#pragma once
/// Node list
typedef vector<uint32_t> Nodelist;
/// Edge structure, neighbor id and the edge weight
typedef pair<uint32_t, double> Edge;
/// Edgelist structure from one source/target node
typedef vector<Edge> Edgelist;
/// Graph structure
typedef vector<Edgelist> Graph;
/// One forward reachable set
typedef vector<size_t> FRset;
/// A set of forward reachable sets
typedef vector<FRset> FRsets;
double decay_factor = 1.0;
/// One reverse reachable set
typedef vector<uint32_t> RRset;
/// A set of reverse reachable sets
typedef vector<RRset> RRsets;
bool optflag;
enum ProbDist
{
    WEIGHTS,
    UNIFORM,
    WC,
    SKEWED,
    PROB_DIST_ERROR
};
enum FuncType
{
    FORMAT,
    IM,
    FUNC_ERROR
};

/// Node element with id and a property value
typedef struct NodeElement
{
    int id;
    double value;
} NodeEleType;

/// Smaller operation for node element
struct smaller
{
    bool operator()(const NodeEleType &Ele1, const NodeEleType &Ele2) const
    {
        return (Ele1.value < Ele2.value);
    }
};

/// Greater operation for node element
struct greater
{
    bool operator()(const NodeEleType &Ele1, const NodeEleType &Ele2) const
    {
        return (Ele1.value > Ele2.value);
    }
};
