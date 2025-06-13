#pragma once
#include <algorithm>
#include <queue>
#include "DataGraph.h"
#include "Define.h"
#define INITIAL -1
// bool compareDegree(const int &a, const int &b);
class TrussDecomposition
{
public:
    // DataGraph *datagraph;
    int kMax;

    map<Edge, int> trussd;
    map<int, set<Edge>> k2edge;
    // map<Edge, int> sup;
    // vector<int> degree;
    // vector<int> total_order;
    // vector<int> order_pointer;
    vector<int> bin;
    TrussDecomposition();
    TrussDecomposition(DataGraph *datagraph);
    void UpdateTrussofInsert2(vector<Edge> &, Edge &e, DataGraph *relSubG);
    // void estimateTruss(Edge &e, DataGraph *relSubG);
    void MaintainSupport_AddEdge(DataGraph *datagraph, vector<set<Edge>> &nonDecSup, Edge e);
    // bool compareDegree(const int, const int);
};

// KTruss::KTruss(DataGraph &dg) : datagraph(dg) {}
// vector<int> KTruss::degree;
// descending sort as trussness
// bool KTruss::compareDegree(const int a, const int b)
// {
//     return degree[a] > degree[b] || (degree[a] == degree[b] && a > b);
// }
