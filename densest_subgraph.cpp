#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <numeric>
#include <unordered_map>
#include <vector>
#include "Define.h"
#include <cstdint>

using namespace std;

void new_edge(int64_t u, int64_t v, double weight, int64_t *to, double *cap, double *flow, int64_t *next, int64_t *fin, long *nEdge)
{
    to[*nEdge] = v;
    cap[*nEdge] = weight;
    flow[*nEdge] = 0;
    next[*nEdge] = fin[u];
    fin[u] = (*nEdge)++;
    to[*nEdge] = u;
    cap[*nEdge] = weight;
    flow[*nEdge] = weight;
    next[*nEdge] = fin[v];
    fin[v] = (*nEdge)++;
}

bool dinic_bfs(int64_t nverts, int64_t src, int64_t dest, int64_t *dist, int64_t *Q, int64_t *fin, int64_t *next, int64_t *to, double *flow, double *cap)
{
    int64_t st, en, i, u, v;
    fill(dist, dist + nverts, -1);
    dist[src] = st = en = 0;
    Q[en++] = src;
    while (st < en)
    {
        u = Q[st++];
        for (i = fin[u]; i >= 0; i = next[i])
        {
            v = to[i];
            if (flow[i] < cap[i] && dist[v] == -1)
            {
                dist[v] = dist[u] + 1;
                Q[en++] = v;
            }
        }
    }
    return dist[dest] != -1;
}

double dinic_dfs(int64_t u, double fl, int64_t src, int64_t dest, int64_t *pro, int64_t *next, int64_t *to, int64_t *dist, double *cap, double *flow)
{
    if (u == dest)
        return fl;
    int64_t v;
    double df;
    for (int64_t &e = pro[u]; e >= 0; e = next[e])
    {
        v = to[e];
        if (flow[e] < cap[e] && dist[v] == dist[u] + 1)
        {
            if (u == src || (cap[e] - flow[e]) <= fl)
            {
                fl = cap[e] - flow[e];
            }
            df = dinic_dfs(v, fl, src, dest, pro, next, to, dist, cap, flow);
            if (df > 0)
            {
                flow[e] += df;
                flow[e ^ 1] -= df;
                return df;
            }
        }
    }
    return 0;
}

void find_cut(int64_t u, int64_t *cut, int64_t *another_pro, int64_t *next, int64_t *to, double *flow, double *cap)
{
    cut[u] = 1;
    for (int64_t &e = another_pro[u]; e >= 0; e = next[e])
    {
        int64_t v = to[e];
        if (flow[e] < cap[e] && cut[v] == 0)
        {
            find_cut(v, cut, another_pro, next, to, flow, cap);
        }
    }
}

double max_flow(double (*edges_info)[3], int64_t nverts, int64_t nedges, int64_t src, int64_t dest, int64_t *Q, int64_t *fin, int64_t *pro, int64_t *dist, int64_t *next, int64_t *to, int64_t *cut, int64_t *another_pro, int64_t *pro3, double *flow, double *cap)
{
    long i;
    fill(fin, fin + nverts, -1);
    fill(cut, cut + nverts, 0);
    long nEdge = 0;
    for (i = 0; i < nedges; i++)
    {
        new_edge(edges_info[i][0], edges_info[i][1], edges_info[i][2], to, cap, flow, next, fin, &nEdge);
    }

    double ret = 0;
    double df;
    while (dinic_bfs(nverts, src, dest, dist, Q, fin, next, to, flow, cap))
    {
        for (i = 0; i < nverts; i++)
        {
            pro[i] = fin[i];
            another_pro[i] = fin[i];
            pro3[i] = fin[i];
        }
        while (true)
        {
            df = dinic_dfs(src, 0, src, dest, pro, next, to, dist, cap, flow);
            if (df)
                ret += df;
            else
                break;
        }
    }
    find_cut(src, cut, another_pro, next, to, flow, cap);
    return ret;
}

void init_edges(double (*edges_info)[3], double *degree, int64_t n, int64_t m, int64_t src, int64_t dest, double g);

/**
 * Compute the densest subgraph given a graph in edge-list format.
 *
 * The edge list format needs to be symmetric, so if (ei[k1],ej[k1]) = (r,s) then there must be (ei[k2],ej[k2]) = (s,r),
 * unless r==s.
 *
 * @param [in] n the number of nodes of the graph
 * @param [in] m the number of edges of the graph (and also the length of arrays ei, ej, w)
 * @param [in] ei a list of sources for each edge (0 <= ei[i] <= n-1) of length m
 * @param [in] ej a list of destinations for each edge (0 <= ei[i] <= n-1) of length m
 * @param [in] w a list of non-negative weights for each edge (0 <= w[i]) of length m
 * @param [out] output a list of vertices in the densest subgraph used, 0 <= output[i] <= n-1 for 0 <= i <= outputlen-1, this array must be capable of holding n
 * @param [in/out] the valid length of the output list and the length of the set of output used.
 * @return the density of the densest subgraph in total_edge_weight/number of vertices
 */
double densest_subgraph(int64_t n, int64_t m, vector<int> &ei, vector<int> &ej, vector<double> &w, vector<int> &output)
{
    int64_t nverts, nedges, src, dest;
    double g, maxflow;
    nverts = n + 2;
    nedges = m + 2 * n;
    int64_t i;
    double(*edges_info)[3] = (double(*)[3])malloc(sizeof(double) * nedges * 3);
    if (edges_info == NULL)
    {
        cout << "malloc fail!" << endl;
    }
    double *degree = (double *)malloc(sizeof(double) * n); // degree of every vertex
    if (degree == NULL)
    {
        cout << "malloc fail!" << endl;
    }
    fill(degree, degree + n, 0);
    unordered_map<int,int> id2seq;
    unordered_map<int,int> seq2id;
    for (i = 0; i < m; i++)
    {
        int u=ei[i],v=ej[i],seq1,seq2;
        if(id2seq.find(u)==id2seq.end()){
            seq1 = id2seq.size();
            id2seq[u]= seq1;
            seq2id[seq1]= u;
        }
        else{
            seq1=id2seq[u];
        }
        if(id2seq.find(v)==id2seq.end()){
            seq2 = id2seq.size();
            id2seq[v]= seq2;
            seq2id[seq2]= v;
        }
        else{
            seq2=id2seq[v];
        }
        edges_info[i][0] = seq1;
        edges_info[i][1] = seq2;
        edges_info[i][2] = w[i];
        degree[(int64_t)edges_info[i][0]] += edges_info[i][2];
        edges_info[i][0]++;
        edges_info[i][1]++;
    }

    double final_degree = 0; // the degree of the densest subgraph
    double L = 0;            // lower bound of Andrew Goldberg's algorithm
    double U = m / 2;        // upper bound of Andrew Goldberg's algorithm
    int64_t *final_cut = (int64_t *)malloc(sizeof(int64_t) * nverts);

    int64_t iter = 0;

    /*malloc enough space to store data used to calculate max flow*/
    int64_t *Q = (int64_t *)malloc(sizeof(int64_t) * nverts);
    int64_t *fin = (int64_t *)malloc(sizeof(int64_t) * nverts);
    int64_t *pro = (int64_t *)malloc(sizeof(int64_t) * nverts);
    int64_t *another_pro = (int64_t *)malloc(sizeof(int64_t) * nverts);
    int64_t *pro3 = (int64_t *)malloc(sizeof(int64_t) * nverts);
    int64_t *dist = (int64_t *)malloc(sizeof(int64_t) * nverts);
    double *flow = (double *)malloc(sizeof(double) * 2 * nedges);
    double *cap = (double *)malloc(sizeof(double) * 2 * nedges);
    int64_t *next = (int64_t *)malloc(sizeof(int64_t) * 2 * nedges);
    int64_t *to = (int64_t *)malloc(sizeof(int64_t) * 2 * nedges);
    int64_t *cut = (int64_t *)malloc(sizeof(int64_t) * nverts);

    /*Andrew Goldberg's algorithm*/
    while (n * (n - 1) * (U - L) >= 1)
    {
        iter++;
        g = (U + L) / 2;
        src = 0;
        dest = nverts - 1;
        init_edges(edges_info, degree, n, m, src, dest, g);
        // cout << "flow iteration " << iter << ": range = (" << L << ", " << U << ")," << "g: " << g << ", solution = ";
        maxflow = max_flow(edges_info, nverts, nedges, src, dest, Q, fin, pro, dist, next, to, cut, another_pro, pro3, flow, cap);
        // cout << maxflow << endl;
        if (accumulate(cut, cut + nverts, 0) == 1)
        {
            U = g;
        }
        else
        {
            L = g;
            for (i = 0; i < nverts; i++)
            {
                final_cut[i] = cut[i];
            }
        }
    }

    final_cut[0] = 0;
    final_cut[nverts - 1] = 0;

    /*retrieve the densest subgraph from the final_cut*/
    int64_t num = 0;
    for (i = 1; i < nverts - 1; i++)
    {
        if (final_cut[i] != 0)
        {
            // output[num++] = i - 1;
            num++;
            output.push_back(seq2id[i - 1]);
            for (int64_t &e = pro3[i]; e >= 0; e = next[e])
            {
                if (final_cut[to[e]] != 0)
                {
                    final_degree += cap[e];
                }
            }
        }
    }
    final_degree /= (2 * num);
    // *outputlen = num;

    /*free space*/
    free(final_cut);
    free(degree);
    degree = NULL;
    free(edges_info);
    free(Q);
    free(fin);
    free(pro);
    free(another_pro);
    free(pro3);
    free(dist);
    free(flow);
    free(cap);
    free(next);
    free(to);
    free(cut);
    return final_degree;
}

void init_edges(double (*edges_info)[3], double *degree, int64_t n, int64_t m, int64_t src, int64_t dest, double g)
{
    long i;
    /*add edge from src to every other vertex*/
    for (i = m; i < m + n; i++)
    {
        edges_info[i][0] = src;
        edges_info[i][1] = i - m + 1;
        edges_info[i][2] = m / 2;
    }
    /*add edge from every other vertex to src*/
    for (i = n + m; i < m + 2 * n; i++)
    {
        edges_info[i][0] = i - m - n + 1;
        edges_info[i][1] = dest;
        edges_info[i][2] = m / 2 + 2 * g - degree[i - m - n];
    }
}
