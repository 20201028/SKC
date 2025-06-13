#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include "Best.h"
#include <stack>
#include <chrono>
#include <algorithm>
#include <math.h>
using namespace std;

void dijkstra(int sourceNode, int targetNode,
              unordered_map<int, unordered_map<int, double>> &graph, vector<int> &path)
{
    // cout << "dijkstra" << endl;
    unordered_map<int, int> predecessor;
    // 使用unordered_map动态记录距离
    unordered_map<int, double> distance;
    unordered_map<int, bool> visited;

    for (auto &e : graph)
    {
        int node = e.first;
        distance[node] = numeric_limits<double>::max();
        visited[node] = false;
    }

    priority_queue<Node, vector<Node>, greater<Node>> pq;

    // 初始化源点
    distance[sourceNode] = 0;
    pq.push({sourceNode, 0});

    while (!pq.empty())
    {
        Node current = pq.top();
        pq.pop();

        int u = current.vertex;

        // 如果当前节点已访问过，跳过
        if (visited[u])
            continue;

        visited[u] = true;

        // 如果到达目标节点，返回最短距离
        if (u == targetNode)
        {
            // return distance[u];
            break;
        }

        // 遍历当前节点的邻居
        if (graph.find(u) != graph.end())
        {
            for (const auto &edge : graph.at(u))
            {
                int v = edge.first;
                double weight = edge.second;

                // 如果找到更短的路径，则更新距离和前驱节点，并将节点加入队列
                if (!visited[v] && distance[u] + weight < distance[v])
                {
                    distance[v] = distance[u] + weight;
                    predecessor[v] = u; // 记录前驱
                    pq.push({v, distance[v]});
                }
            }
        }
    }

    // 如果目标节点不可达，返回 -1 表示无解
    // return -1;
    int current = targetNode;

    // 回溯前驱数组，构建路径
    while (predecessor.find(current) != predecessor.end())
    {
        if (current != targetNode)
            path.push_back(current);
        current = predecessor.at(current);
    }
}

void GDS(double &maxRho, set<int> &maxRV, set<int> &maxR, unordered_map<int, double> &vert2weightDs,
         unordered_map<int, double> &vert2weightCC, int &gCCId, set<int> &RCCSubGDS,
         unordered_map<int, unordered_map<int, double>> &simG, DataGraph *relSubG,
         map<int, vector<Edge>, greater<int>>::iterator rel, unordered_map<int, set<int>> &rel2verts)
{
    int r = rel->first;

    unordered_map<int, set<int>> ccId2RelVerts;

    auto start = std::chrono::high_resolution_clock::now();
    for (auto cc : relSubG->CC)
    {
        // 与每个连通分量交
        set<int> verts;
        for (int v : cc.second)
        {
            verts.insert(relSubG->seq2id[v]);
        }
        set<int> comVerts;
        set_intersection(rel2verts[r].begin(), rel2verts[r].end(), verts.begin(), verts.end(), inserter(comVerts, comVerts.begin()));
        if (!comVerts.empty())
        {
            ccId2RelVerts.emplace(cc.first, comVerts);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // cout << "Found the connected components containing the relation. " << duration.count() << " milliseconds." << endl;
    // connected components

    int i = 0;
    int canNum = 0;
    // for each connected component, find a cand, and return the max rho, max rho vertices, and simG
    // double maxRho = 0;
    for (auto cc : ccId2RelVerts)
    {
        // if(i==29){
        //     bool flag = false;
        // }
        // cout << "Connected component " << i++ << ": " << endl;
        int ccId = cc.first;
        set<int> curCC = relSubG->CC[ccId];
        set<int> R;
        for (auto &v : cc.second)
            R.insert(relSubG->id2seq[v]);
        /////////////////////////density of whole graph
        int vN = curCC.size();
        unordered_map<int, double> vert2weight;
        set<int> removeV;

        // record weight for each vertices in this connected component;
        // double maxWeight = 0;
        start = std::chrono::high_resolution_clock::now();
        for (auto &v : curCC)
        {
            double w = 0;
            auto it = simG.find(v);
            if (it != simG.end())
                for (auto &u : simG[v])
                {
                    if (!curCC.count(u.first))
                        continue;
                    w += u.second + denStyle;
                }
            vert2weight.emplace(v, w);

            // if(w <= maxRho) lessMRhoV.push(v);
            // maxWeight = max(maxWeight, w);
        }
        unordered_map<int, double> tempV2WCC = vert2weight;
        set<int> tempRCCSubGDS = R;
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // cout << "maxWeight:" << maxWeight << ", " << duration.count() << " milliseconds." << endl;
        // unordered_map<int, double> vert2weightCopy = vert2weight;
        // if (maxWeight <= 2 * maxRho)
        // {
        //     cout << maxWeight / 2 << " <= " << maxRho << endl;
        //     cout << "first pruning is successful" << endl;
        //     continue;
        // }

        double rho = 0;
        for (auto &v : vert2weight)
            rho += v.second;
        rho /= 2 * vN;
        int z = -1;

        double tempMaxRho = rho;
        set<int> tempMaxVs = curCC;
        set<int> tempMaxR = R;

        unordered_map<int, double> tempVert2weight = vert2weight; // in there, delete vertices's weight, but in basic do not delete vertices, just do not update those vertices not in result graph
        // cout << "Rho: " << tempMaxRho << endl;

        while (curCC.size() > 1)
        {
            auto it = curCC.begin();
            int vMin1 = *it++, vMin2 = *it++; // the minimum value is vMin1
            if (vert2weight[vMin1] > vert2weight[vMin2])
            {
                swap(vMin1, vMin2);
            }

            for (; it != curCC.end(); it++)
            {
                int v = *it;
                double vW = vert2weight[v];
                if (vW < vert2weight[vMin1])
                {
                    vMin2 = vMin1;
                    vMin1 = v;
                }
                else if (vW < vert2weight[vMin2])
                {
                    vMin2 = v;
                }
            }
            int v = vMin1;
            if (R.size() == 1 && R.count(v))
            {
                z = v;
                curCC.erase(v);
                v = vMin2;
            }

            removeV.insert(v);
            curCC.erase(v);

            // update weight
            // total += vert2weight[v];
            for (auto &e : simG[v])
            {
                if (!curCC.count(e.first) && z != e.first)
                    continue;
                // simG[e.first].erase(v);
                double t = vert2weight[e.first];

                vert2weight[e.first] = t - (e.second + denStyle);
            }
            // vert2weight.erase(v);
            // simG.erase(v);
            rho = (rho * vN - vert2weight[v]) / --vN;
            vert2weight.erase(v);
            R.erase(v);

            if (rho > tempMaxRho)
            {
                // cout << "rho: " << rho << endl;
                tempMaxRho = rho;
                tempMaxVs = curCC;
                if (z != -1 && !curCC.count(z))
                    tempMaxVs.insert(z);
                tempMaxR = R;
                tempVert2weight = vert2weight;
            }
        }

        // return the maximum from each cc.
        if (tempMaxRho > maxRho)
        {
            // cout << "Found a candidate, rho = " << tempMaxRho << endl;
            maxRho = tempMaxRho;
            maxRV = tempMaxVs;
            maxR = tempMaxR;
            gCCId = ccId;
            vert2weightDs = tempVert2weight;
            vert2weightCC = tempV2WCC;
            RCCSubGDS = tempRCCSubGDS;
        }
    }
}
void connect(int node, unordered_set<int> &visited, unordered_map<int, unordered_map<int, int>> &G)
{
    visited.insert(node);

    for (auto &e : G[node])
    {
        int neighbor = e.first;
        if (visited.find(neighbor) == visited.end())
        {
            connect(neighbor, visited, G);
        }
    }
}
void GKMDS(int &kMax, set<int> &maxRhoG, double &maxRho, set<int> R,
           unordered_map<int, double> vert2weight, unordered_map<int, unordered_map<int, double>> &simG,
           unordered_map<int, unordered_map<int, int>> edge2sup)
{
    // get the support of each edge in maxRhoG
    // unordered_map<int, unordered_map<int, int>> edge2sup;
    set<int> kMaxV;
    map<int, vector<Edge>> sup2edge;
    // for (auto U : edge2sup)
    // {
    //     int u = U.first;
    //     for (auto V : U.second)
    //     {
    //         int v = V.first;
    //         if (v > u)
    //             continue;
    //         Edge e = make_pair(v, u);
    //         int sup = V.second;
    //         if (!sup2edge.count(sup))
    //         {
    //             sup2edge.emplace(sup, vector<Edge>{e});
    //         }
    //         else
    //             sup2edge[sup].push_back(e);
    //     }
    // }
    for (auto u : maxRhoG)
    {
        for (auto V : edge2sup[u])
        {
            int v = V.first;
            if (v > u)
                continue;
            Edge e = make_pair(v, u);
            int sup = V.second;
            if (!sup2edge.count(sup))
            {
                sup2edge.emplace(sup, vector<Edge>{e});
            }
            else
                sup2edge[sup].push_back(e);
        }
    }

    while (!sup2edge.empty())
    {
        // relSubG->support(edge2sup, sup2edge, maxRhoG);
        kMaxV = maxRhoG; // store the kMax vertex set
        int vN = maxRhoG.size();
        auto minSupE = sup2edge.begin();
        kMax = minSupE->first + 2; // the initial kMaxV is the maxRhoG
        // cout << "first kMax: " << kMax << endl;
        // set<Edge> visit;
        unordered_set<int> visit;
        unordered_map<int, set<int>> neis;
        while (!minSupE->second.empty())
        {
            Edge e = minSupE->second.back();
            minSupE->second.pop_back();
            int u = e.first;
            int v = e.second;
            // if (!edge2sup.count(u) || !edge2sup[u].count(v))
            //     continue;

            edge2sup[u].erase(v); // remove edge
            edge2sup[v].erase(u);
            visit.insert(u);
            visit.insert(v);
            for (auto w2sup : edge2sup[u])
            {
                int w = w2sup.first;
                if (edge2sup[v].count(w))
                {
                    Edge e1 = make_pair(min(w, u), max(w, u));
                    int oldsup = edge2sup[u][w];
                    int uwSup = --edge2sup[u][w];
                    --edge2sup[w][u];
                    // vector<Edge> &uw2sup = sup2edge[uwSup+1];
                    // auto ituw = find(uw2sup.begin(),uw2sup.end(), e1);
                    // uw2sup.erase(ituw);

                    if (uwSup == kMax - 2) // if < , which means the oldsup=14, it already in minsupE
                    {
                        minSupE->second.push_back(e1);
                    }
                    else if (uwSup > kMax - 2)
                    {
                        auto newSup = sup2edge.find(uwSup);
                        if (newSup == sup2edge.end())
                        {
                            sup2edge.emplace(uwSup, vector<Edge>{e1});
                        }
                        else
                            newSup->second.push_back(e1);
                    }
                    if (uwSup >= kMax - 2)
                    {
                        auto it1 = find(sup2edge[oldsup].begin(), sup2edge[oldsup].end(), e1);
                        sup2edge[oldsup].erase(it1);
                    }

                    Edge e2 = make_pair(min(w, v), max(w, v));
                    oldsup = edge2sup[v][w];
                    int vwSup = --edge2sup[v][w];
                    --edge2sup[w][v];
                    if (vwSup == kMax - 2)
                    {
                        minSupE->second.push_back(e2);
                    }
                    else if (vwSup > kMax - 2)
                    {
                        auto newSup = sup2edge.find(vwSup);
                        if (newSup == sup2edge.end())
                        {
                            sup2edge.emplace(vwSup, vector<Edge>{e2});
                        }
                        else
                            newSup->second.push_back(e2);
                    }
                    if (vwSup >= kMax - 2)
                    {
                        auto it2 = find(sup2edge[oldsup].begin(), sup2edge[oldsup].end(), e2);
                        sup2edge[oldsup].erase(it2);
                    }
                }
            }
        }

        // when there exists a larger k
        sup2edge.erase(kMax - 2); // sup2edge's edges corresponding to each sup only add, not delete

        set<int> removeV;
        for (int u : visit)
        {
            if (edge2sup[u].empty())
            {
                removeV.insert(u);
                edge2sup.erase(u); // delete isolated vertices
            }
        }
        if (edge2sup.empty())
            return;
        // sort(removeV.begin(), removeV.end());
        set<int> comV;
        set_difference(R.begin(), R.end(), removeV.begin(), removeV.end(), inserter(comV, comV.begin()));
        if (comV.empty())
            return; // means that there exists a vertex can be removed because at least have a vertex's relevance is r. Thus, the subgraph's trussness is lesser than and equal before.

        // check the connective after deleting vertices
        unordered_set<int> visited;
        int startV = edge2sup.begin()->first;
        connect(startV, visited, edge2sup);
        if (visited.size() != edge2sup.size())
            return;
        vector<int> lERho;
        for (int v : removeV)
        {
            if (vert2weight[v] <= maxRho)
                lERho.push_back(v);
        }
        visited.clear();
        double total = 0;
        while (!lERho.empty())
        {
            int v = lERho.back();
            visited.insert(v);
            lERho.pop_back();
            for (auto &e : simG[v])
            {
                if (!maxRhoG.count(e.first) || visited.count(e.first))
                    continue;
                // simG[e.first].erase(v);
                total += (e.second + denStyle);
                double t = vert2weight[e.first];
                vert2weight[e.first] = t - (e.second + denStyle);
                if (vert2weight[e.first] <= maxRho && removeV.count(e.first))
                {
                    lERho.push_back(e.first);
                }
            }
        }
        if (visited.size() != removeV.size())
            return;
        // next loop
        for (int v : removeV)
        {
            R.erase(v);
            maxRhoG.erase(v);
        }
        // recompute the rho
        maxRho = (maxRho * vN - total) / (maxRhoG.size());
    }
    maxRhoG.clear();
    maxRhoG = kMaxV;
}

bool compareByDensity(const std::pair<int, double> &a, const std::pair<int, double> &b)
{
    return a.second > b.second;
}
bool constructGraph(int &src, int &dst, int gCCId, DataGraph *G, unordered_map<int, double> &vert2weight, unordered_map<int, unordered_map<int, double>> &g, set<int> &a, set<int> &b)
{
    set<int> processedNodes; // 记录已经处理的点
    // 将 a 集合中的点合并为一个新点
    // int a_node = -1; // 用 -1 表示 a 的合并点
    for (int u : a)
    {

        for (int v : G->AdjList[u]) // 因为是对G的一个连通分量进行构图，所以不需要判断v是否在G->CC[gCCId]中
        {
            if (a.count(v))
                continue; // 忽略 a 内部的连接

            if (processedNodes.count(v))
                continue;
            if (b.count(v))
                return true;
            processedNodes.insert(v); // 记录已处理点
            g[v][src] = 0;
            g[src][v] = exp(-vert2weight.at(v) - G->AdjList[v].size());
        }
    }

    // 将 b 集合中的点合并为一个新点
    // int b_node = -2; // 用 -2 表示 b 的合并点
    processedNodes.clear();
    for (int u : b)
    {
        // processedNodes.insert(u); // 记录已处理点
        for (int v : G->AdjList[u])
        {
            if (b.count(v))
                continue; // 忽略 b 内部的连接
            if (processedNodes.count(v))
                continue;
            processedNodes.insert(v);
            g[v][dst] = 0;
            g[dst][v] = exp(-vert2weight.at(v) - G->AdjList[v].size());
        }
    }

    // 保留其余点和边，转换为有向图
    for (int u : G->CC[gCCId])
    {
        if (a.count(u) || b.count(u))
            continue; // 跳过已处理点
        for (int v : G->AdjList[u])
        {
            if (a.count(v) || b.count(v))
                continue;
            g[v][u] = exp(-vert2weight.at(u) - G->AdjList[u].size()); // 使用 u 的权重作为边权
        }
    }
    return false;
}
void expand(SkyGroupCand &c, double &expRho, set<int> &expRV, set<int> &expR, unordered_map<int, unordered_map<int, double>> &simG,
            unordered_map<int, unordered_map<int, int>> &expSup, DataGraph *subG, set<int> &RCCSubGDS,
            unordered_map<int, double> &V2WCCSubGDS, unordered_map<int, double> &expW, int k, int &kMax)
{

    // cout << "expand" << endl;
    // maxVert2weight = aW;
    // get vertices whose incident edges' support less than k-2
    auto start = std::chrono::high_resolution_clock::now();
    unordered_map<int, int> gap;
    for (auto &U : expSup)
    {
        int u = U.first;
        int minSup = U.second.begin()->second;
        for (auto &V : U.second)
        {
            int v = V.first;
            int support = V.second;
            minSup = min(minSup, support);
        }
        if (minSup < k - 2)
        {
            gap.emplace(u, minSup);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "gap : " << duration.count() << " milliseconds." << std::endl;
    // get vertices connected to the vertices whose minimum incident edges' support less than k-2
    unordered_map<int, double> nei;
    // map<int, double> nei;
    // for (auto &U : gap)
    // {
    //     int u = U.first;
    //     for (int v : subG->AdjList[u])
    //     {
    //         if (!expRV.count(v))
    //         {
    //             auto iter = nei.find(v);
    //             if (iter != nei.end())
    //                 iter->second++;
    //             else
    //             {
    //                 nei.emplace(v, 1);
    //             }
    //         }
    //     }
    // }
    // gap + weighted 1121033>gap&weighted 737445>gap 710588>gap + nei( + weighted)
    start = std::chrono::high_resolution_clock::now();
    for (auto &U : gap)
    {
        int u = U.first;
        for (int v : subG->AdjList[u])
        {
            if (!expRV.count(v))
            {
                auto iter = nei.find(v);
                if (iter != nei.end())
                    iter->second++;
                else
                {
                    nei.emplace(v, V2WCCSubGDS[v]);
                    // nei.emplace(v, 1);
                }
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "nei : " << duration.count() << " milliseconds." << std::endl;
    // gap + nei + weighted in expanded subgraph: if considering nei, the speed to k is very slow
    // for (int u : expRV)
    // {
    //     for (int v : subG->AdjList[u])
    //     {
    //         if (!expRV.count(v))
    //         {
    //             auto iter = nei.find(v);
    //             if (iter != nei.end())
    //                 iter->second++;
    //             else
    //             {
    //                 nei.emplace(v, V2WCCSubGDS[v]);
    //                 // nei.emplace(v, 1);
    //             }
    //             if (gap.count(u))
    //                 nei[v]++;
    //         }
    //     }
    // }
    // selected vertices who maximum the sum between the number of gap neighbor and maximum weight
    if (nei.empty())
        return;
    int maxScoreV = nei.begin()->first;
    double maxScore = nei.begin()->second;
    for (auto &vscore : nei)
    {
        if (vscore.second > maxScore)
        {
            maxScoreV = vscore.first;
            maxScore = vscore.second;
        }
    }
    // // 将 unordered_map 转换为 vector
    // std::vector<std::pair<int, double>> vec(nei.begin(), nei.end());

    // // 自定义比较函数，按值降序排序
    // auto compare = [](const std::pair<int, double> &a, const std::pair<int, double> &b)
    // {
    //     return a.second > b.second;
    // };

    // // 对 vector 进行排序
    // std::sort(vec.begin(), vec.end(), compare);
    // int maxScoreV = vec.begin()->first;
    // double maxScore = vec.begin()->second;
    // double maxW = V2WCCSubGDS[maxScoreV]+subG->AdjList[maxScoreV].size();
    // for (const auto &pair : vec)
    // {
    //     if (pair.second != maxScore)
    //         break;
    //     if (V2WCCSubGDS[pair.first] > maxW)
    //     {
    //         maxScoreV = pair.first;
    //         // maxW = V2WCCSubGDS[pair.first]+subG->AdjList[pair.first].size();
    //         maxW = V2WCCSubGDS[pair.first];
    //         maxW = V2WCCSubGDS[pair.first]+subG->AdjList[pair.first].size();
    //     }
    // }

    expSup.emplace(maxScoreV, unordered_map<int, int>());
    start = std::chrono::high_resolution_clock::now();
    unordered_set<int> vNeiExpG;
    for (int &u : subG->AdjList[maxScoreV])
    {
        if (expRV.count(u))
        {
            vNeiExpG.insert(u);
        }
    }
    // update support
    unordered_set<int> visit;
    for (int u : vNeiExpG)
    {
        visit.insert(u);

        int sup = 0;
        for (auto &W : expSup[u])
        {
            int w = W.first;
            if (vNeiExpG.count(w))
            {
                sup++;
                if (!visit.count(w))
                {
                    expSup[u][w]++;
                    expSup[w][u]++;
                }
            }
        }
        expSup[maxScoreV][u] = sup;
        expSup[u][maxScoreV] = sup;
    }
    // updata weight
    double tw = 0;
    for (auto &U : simG[maxScoreV])
    {
        int u = U.first;

        if (expRV.count(u))
        {
            expW[u] += U.second;
            tw += U.second;
        }
    }
    expW.emplace(maxScoreV, tw);
    // update vertices
    expRV.insert(maxScoreV); // cannot change the origin a
    // vector<int> vertices(expRV.begin(), expRV.end());
    // DataGraph *sub = subG->getSubG(vertices);
    // TrussDecomposition *td = new TrussDecomposition(sub);
    // for (auto &kt : td->k2edge)
    // {
    //     if (kt.second.size() != 0)
    //     {
    //         cout << "k value: " << kt.first << endl;
    //         break;
    //     }
    // }

    // update rho
    expRho = (expRho * (expRV.size() - 1) + tw) / expRV.size();
    // update R
    if (RCCSubGDS.count(maxScoreV))
        expR.insert(maxScoreV);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "update : " << duration.count() << " milliseconds." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    double maxRho = expRho;
    set<int> maxRV = expRV;
    // set<int> maxR = expR;
    unordered_map<int, unordered_map<int, int>> maxSup = expSup;
    // int kMax = 0;
    GKMDS(kMax, maxRV, maxRho, expR, expW, simG, maxSup);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "GKMDS : " << duration.count() << " milliseconds." << std::endl;
    if (kMax >= k)
    {
        c.k = kMax;
        c.rho = maxRho;
        for (int v : maxRV)
        {
            c.vertices.insert(subG->seq2id[v]);
        }

        return;
    }
    else
        expand(c, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, k, kMax);
}
void Best(map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
    auto rel = rel2edge.begin();
    DataGraph *relSubG = new DataGraph();
    vector<SkyGroupCand> cand;

    ///////////////////////////////////////
    while (rel != rel2edge.end())
    {
        auto start = std::chrono::high_resolution_clock::now();
        cand.clear();
        for (auto &e : rel->second)
        {
            int src = e.first;
            int dst = e.second;
            // relSubG->addEdgeAndMainConnect(src, dst, simE, simG);
            relSubG->addEdgeNoMatinC(src, dst);
        }

        // if (rel->first == 2)
        // {
        //     ++rel;
        //     continue;
        // }
        // int curK = 3;

        TrussDecomposition *kt = new TrussDecomposition(relSubG);
        // cout << "TrussMaintance completed" << endl;
        int curK = 1;
        // int curK = kt->kMax-1;
        // find k+1-truss
        while (curK < kt->kMax)
        {

            DataGraph *subG = new DataGraph();
            unordered_map<int, unordered_map<int, double>> simG;
            for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
            {
                int k = it->first;
                if (k <= curK)
                {
                    it = kt->k2edge.erase(it);
                }
                else
                {
                    for (auto &e : it->second)
                    {
                        int src = relSubG->seq2id[e.first];
                        int dst = relSubG->seq2id[e.second];
                        subG->addEdgeAndMainConnect(src, dst, simE, simG);
                    }

                    ++it;
                }
            }
            double rhoDS = 0;
            set<int> vDS;
            set<int> RDS;
            unordered_map<int, double> v2wDS;
            unordered_map<int, double> V2WCCSubGDS;
            set<int> RCCSubGDS;
            int ccIdSubG = 0;
            // find the densest subgraph, get the maxRho, maxRV, maxR, maxVert2weight, gCCId
            GDS(rhoDS, vDS, RDS, v2wDS, V2WCCSubGDS, ccIdSubG, RCCSubGDS, simG, subG, rel, rel2verts);
            if (vDS.empty())
                break;
            // gccid record the id of the connected component from densest subgraph
            // create a new subgraph
            unordered_map<int, set<int>> ccDS;
            // DataGraph *sub = new DataGraph();
            // if (rel->first == 1 && curK == 3)
            // {
            //     bool flag = true;
            // }
            subG->getCC(ccDS, vDS);
            // subG->getSubG(ccDS, DS, maxRV);
            int kMax = 0;
            // if (curK == 15)
            // {
            //     bool flag = true;
            // }
            if (ccDS.size() == 1)
            {
                // cout << "Only one connected component" << endl;
                unordered_map<int, unordered_map<int, int>> DS; // the densest subgraph
                subG->support(DS, vDS);
                // vector<int> vertices(vDS.begin(), vDS.end());
                // DataGraph *sub = subG->getSubG(vertices);
                // TrussDecomposition *td = new TrussDecomposition(sub);
                double maxRho = rhoDS;
                set<int> maxRV = vDS;
                // set<int> maxR = RDS;
                // unordered_map<int, double> maxVert2weight = v2wDS;
                // get kMax-truss holding maxRho or larger rho
                GKMDS(kMax, maxRV, maxRho, RDS, v2wDS, simG, DS); // RDS, DS and v2wDS不需要返回新值，所以可以不单独赋值，传入函数的是副本
                SkyGroupCand c = SkyGroupCand();
                if (kMax < curK + 1)
                {
                    // cout << "curK: " << curK << endl;
                    expand(c, rhoDS, vDS, RDS, simG, DS, subG, RCCSubGDS, V2WCCSubGDS, v2wDS, curK + 1, kMax);
                    // output kMax, maxRho, maxRV
                    c.rel = rel->first;
                    // cand.push_back(c);
                }
                else
                {
                    // is result, output kMax, maxRho, maxRV

                    c.rel = rel->first;
                    c.k = kMax;
                    c.rho = maxRho;
                    for (int v : maxRV)
                    {
                        c.vertices.insert(subG->seq2id[v]);
                    }
                }
                curK = c.k;
                cand.push_back(c);
            }
            else
            {
                // cout << "Some connected components" << endl;
                // 连通分量可能不包含rel点，所以从包含rel的分量中找到rho最大的分量，判断最大k值是否大于k，如果不是则拓展。然后从剩余的分量里找最大，权重最短路径连接。
                // 用2-近似得到密度最大，而不用k-truss密度最大是因为k-truss不是单调的，通过剥蒜式缩减得到的密度可能不是最大的，即其子图中可能还包含一个密度更大的k-truss。
                // unordered_map<int, double> ccDs2Rho;
                vector<pair<int, double>> ccDs2Rho;
                unordered_map<int, unordered_map<int, double>> ccDs2vert2weight;
                bool isSingle = true;
                for (auto c : ccDS)
                {
                    if (c.second.size() > 1)
                        isSingle = false;
                }
                if (isSingle)
                {
                    for (auto c : ccDS)
                    {
                        ccDs2vert2weight.emplace(c.first, unordered_map<int, double>{{*c.second.begin(), 0}});
                        ccDs2Rho.push_back(make_pair(c.first, 0));
                    }
                    sort(ccDs2Rho.begin(), ccDs2Rho.end(), [&v2wDS](const pair<int, double> &a, const pair<int, double> &b)
                         { return v2wDS[a.first] > v2wDS[b.first]; });
                }
                else
                {
                    for (auto c : ccDS)
                    {
                        ccDs2vert2weight.emplace(c.first, unordered_map<int, double>());
                        auto &vert2weight = ccDs2vert2weight[c.first];
                        // whether contain rel vertices
                        set<int> &curCC = c.second;
                        double total = 0;
                        for (int v : curCC)
                        {
                            double w = 0;
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : simG[v])
                                {
                                    if (!curCC.count(u.first))
                                        continue;
                                    w += u.second + denStyle;
                                }
                            vert2weight.emplace(v, w);
                            total += w;
                        }
                        total /= (2 * curCC.size());
                        ccDs2Rho.push_back(make_pair(c.first, total));
                    }
                    sort(ccDs2Rho.begin(), ccDs2Rho.end(), compareByDensity);
                }

                // find the max cc containing rel
                double aRho = 0;
                set<int> RA;
                int ccRDsmaxRho;
                for (auto &c : ccDs2Rho)
                {
                    int ccId = c.first;
                    set<int> &curCC = ccDS[ccId];
                    set<int> com;
                    set_intersection(curCC.begin(), curCC.end(), RDS.begin(), RDS.end(), inserter(com, com.end()));

                    if (!com.empty())
                    {
                        ccRDsmaxRho = ccId;
                        RA = com;

                        aRho = c.second;
                        break;
                    }
                }
                set<int> &a = ccDS[ccRDsmaxRho];
                unordered_map<int, double> &aW = ccDs2vert2weight[ccRDsmaxRho];
                double maxRho = aRho;
                set<int> maxRV = a;
                // set<int> maxR = ccRDs;
                unordered_map<int, unordered_map<int, int>> aSup;
                SkyGroupCand tempc = SkyGroupCand(); // store the current densest k-truss
                if (!isSingle)
                {

                    // check wether is k-truss if not compute the gap between trussness of vertices in the max cc and the expected k
                    subG->support(aSup, a); // only the edges in cc of maxRV have been counted
                    // 判断是否满足k约束，这里传入的是副本，不对原变量进行修改，如果满足，返回k值、子图、密度
                    // vector<int> vertices(a.begin(), a.end());
                    // DataGraph *sub = subG->getSubG(vertices);
                    // TrussDecomposition *td = new TrussDecomposition(sub);
                    // for (auto &kt : td->k2edge)
                    // {
                    //     if (kt.second.size() != 0)
                    //     {
                    //         cout << "k value: " << kt.first << endl;
                    //         break;
                    //     }
                    // }
                    GKMDS(kMax, maxRV, maxRho, RA, aW, simG, aSup);
                    // SkyGroupCand tempc = SkyGroupCand(); // store the current densest k-truss
                    if (kMax < curK + 1)
                    {
                        // double expRho = aRho;
                        // set<int> expRV = a;
                        // set<int> expR = RA;
                        // unordered_map<int, unordered_map<int, int>> expSup = aSup;
                        // unordered_map<int, double> expW = aW;
                        // expand(tempc, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, curK + 1, kMax); // cannot change a
                        expand(tempc, aRho, a, RA, simG, aSup, subG, RCCSubGDS, V2WCCSubGDS, aW, curK + 1, kMax);
                    }
                    else
                    {
                        // is a candidate
                        // c.rel = rel->first;
                        tempc.k = kMax;
                        tempc.rho = maxRho;
                        for (int v : maxRV)
                        {
                            tempc.vertices.insert(subG->seq2id[v]);
                        }
                        // cand.push_back(c);
                    }
                }
                for (auto &c : ccDs2Rho)
                {
                    if (c.first == ccRDsmaxRho)
                        continue;
                    // use least high weight vertices in cc of subG to connect the current cc with the max cc in the rest of ccs, and loop again.
                    // first check whether the expend a connects with b, if so, do not need to connect again.
                    set<int> &b = ccDS[c.first];
                    set<int> bRest;
                    set_difference(b.begin(), b.end(), a.begin(), a.end(), inserter(bRest, bRest.end()));
                    if (bRest.empty())
                        continue;
                    else if (bRest.size() == b.size())
                    {
                        unordered_map<int, unordered_map<int, double>> transG;
                        vector<int> path;
                        int src = -1, dst = -2;

                        bool connt = constructGraph(src, dst, ccIdSubG, subG, V2WCCSubGDS, transG, a, b); // construct graph on the subgraph containing DS
                        if (!connt)
                            dijkstra(src, dst, transG, path); // get added vertices
                        // cout << "path: " << path.size() << endl;
                        // calculate weight of new connected component
                        unordered_map<int, double> &bW = ccDs2vert2weight[c.first];
                        for (auto &vw : bW)
                        {
                            aW.emplace(vw);
                        }
                        aRho = aRho * a.size() + c.second * b.size();
                        for (int v : b)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!a.count(u.first))
                                        continue;
                                    aW[u.first] += u.second + denStyle;
                                    aW[v] += u.second + denStyle;
                                    aRho += u.second + denStyle;
                                }
                        }
                        a.insert(b.begin(), b.end());

                        double tt = 0;
                        // weight sum of shortest path
                        for (int v : path)
                        {
                            double w = 0;
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!count(path.begin(), path.end(), u.first))
                                        continue;
                                    w += u.second + denStyle;
                                }
                            aW.emplace(v, w);
                            tt += w;
                        }
                        aRho += tt / 2;
                        for (int v : path)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!a.count(u.first))
                                        continue;
                                    aW[u.first] += u.second + denStyle;
                                    aW[v] += u.second + denStyle;
                                    aRho += u.second + denStyle;
                                }
                        }
                        for (int v : path)
                        {
                            a.insert(v);
                        }
                        aRho /= a.size();

                        // update R
                        for (int v : RCCSubGDS)
                        {
                            if (count(path.begin(), path.end(), v))
                                RA.insert(v);
                            if (b.count(v))
                                RA.insert(v);
                        }
                        b.clear();
                    }
                    else // already connected
                    {

                        for (int v : bRest)
                        {
                            aW.emplace(v, 0);
                        }
                        aRho = aRho * a.size();
                        for (int v : bRest)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &U : it->second)
                                {
                                    int u = U.first;
                                    if (a.count(u))
                                    {
                                        aW[u] += U.second + denStyle;
                                        aW[v] += U.second + denStyle;
                                        aRho += U.second + denStyle;
                                    }
                                    else if (v < u && bRest.count(u))
                                    {
                                        aW[v] += U.second + denStyle;
                                        aW[u] += U.second + denStyle;
                                        aRho += U.second + denStyle;
                                    }
                                }
                        }
                        a.insert(bRest.begin(), bRest.end());
                        aRho /= a.size();

                        // update R
                        for (int v : RCCSubGDS)
                        {
                            if (b.count(v))
                                RA.insert(v);
                        }
                    }
                    maxRho = aRho;
                    maxRV = a;
                    // maxVert2weight = aW;
                    // update support
                    aSup.clear();
                    subG->support(aSup, a); // only the edges in cc of maxRV have been counted
                    // vector<int> vertices(a.begin(), a.end());
                    // DataGraph *sub = subG->getSubG(vertices);
                    // TrussDecomposition *td = new TrussDecomposition(sub);
                    // for (auto &kt : td->k2edge)
                    // {
                    //     if (kt.second.size() != 0)
                    //     {
                    //         cout << "k value: " << kt.first << endl;
                    //         break;
                    //     }
                    // }
                    GKMDS(kMax, maxRV, maxRho, RA, aW, simG, aSup);
                    // expend the max cc.
                    if (kMax < curK + 1)
                    {
                        // double expRho = aRho;
                        // set<int> expRV = a;
                        // set<int> expR = RA;
                        // unordered_map<int, unordered_map<int, int>> expSup = aSup;
                        // unordered_map<int, double> expW = aW;
                        // SkyGroupCand expCand = SkyGroupCand();
                        // expand(expCand, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, curK + 1, kMax);
                        // if (expCand.rho >= tempc.rho)
                        // {
                        //     tempc.k = expCand.k;
                        //     tempc.rho = expCand.rho;
                        //     tempc.vertices.clear();
                        //     for (int v : expCand.vertices)
                        //     {
                        //         tempc.vertices.insert(subG->seq2id[v]);
                        //     }
                        // }
                        SkyGroupCand expCand = SkyGroupCand();
                        expand(expCand, aRho, a, RA, simG, aSup, subG, RCCSubGDS, V2WCCSubGDS, aW, curK + 1, kMax);
                        if (expCand.rho >= tempc.rho)
                        {
                            tempc.k = expCand.k;
                            tempc.rho = expCand.rho;
                            tempc.vertices = expCand.vertices;
                            // tempc.vertices.clear();
                            // for (int v : expCand.vertices)
                            // {
                            //     tempc.vertices.insert(subG->seq2id[v]);
                            // }
                        }
                    }
                    else
                    {
                        if (maxRho >= tempc.rho)
                        {
                            tempc.k = kMax;
                            tempc.rho = maxRho;
                            tempc.vertices.clear();
                            for (int v : maxRV)
                            {
                                tempc.vertices.insert(subG->seq2id[v]);
                            }
                        }
                    }
                }
                tempc.rel = rel->first;
                cand.push_back(tempc);
                // cout << "k: " << curK << ", K: " << tempc.k << endl;
                curK = tempc.k;
            }
        }

        // cout << "Find candidate skyline groups" << endl;

        skyline(result, cand);
        // for (const auto &point : cand)
        // {
        //     cout << "rel: " << point.rel << " k: " << point.k << " rho: " << point.rho << " maxRV: ";
        //     for (auto v : point.vertices)
        //         cout << v << ", ";
        //     cout << endl;
        // }

        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "skyline2 : " << duration.count() << " milliseconds." << std::endl;
        ++rel;
        delete kt;
    }

    delete relSubG;
}
void Best(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
    auto rel = rel2edge.begin();
    DataGraph *relSubG = new DataGraph();
    vector<SkyGroupCand> cand;

    ///////////////////////////////////////
    while (rel != rel2edge.end())
    {
        auto start = std::chrono::high_resolution_clock::now();
        cand.clear();
        for (auto &e : rel->second)
        {
            int src = e.first;
            int dst = e.second;
            // relSubG->addEdgeAndMainConnect(src, dst, simE, simG);
            relSubG->addEdgeNoMatinC(src, dst);
        }

        // if (rel->first == 2)
        // {
        //     ++rel;
        //     continue;
        // }
        // int curK = 3;

        TrussDecomposition *kt = new TrussDecomposition(relSubG);
        // cout << "TrussMaintance completed" << endl;
        int curK = k - 1;
        // int curK = kt->kMax-1;
        // find k+1-truss
        while (curK < kt->kMax)
        {

            DataGraph *subG = new DataGraph();
            unordered_map<int, unordered_map<int, double>> simG;
            for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
            {
                int k = it->first;
                if (k <= curK)
                {
                    it = kt->k2edge.erase(it);
                }
                else
                {
                    for (auto &e : it->second)
                    {
                        int src = relSubG->seq2id[e.first];
                        int dst = relSubG->seq2id[e.second];
                        subG->addEdgeAndMainConnect(src, dst, simE, simG);
                    }

                    ++it;
                }
            }
            double rhoDS = 0;
            set<int> vDS;
            set<int> RDS;
            unordered_map<int, double> v2wDS;
            unordered_map<int, double> V2WCCSubGDS;
            set<int> RCCSubGDS;
            int ccIdSubG = 0;
            // find the densest subgraph, get the maxRho, maxRV, maxR, maxVert2weight, gCCId
            GDS(rhoDS, vDS, RDS, v2wDS, V2WCCSubGDS, ccIdSubG, RCCSubGDS, simG, subG, rel, rel2verts);
            if (vDS.empty())
                break;
            // gccid record the id of the connected component from densest subgraph
            // create a new subgraph
            unordered_map<int, set<int>> ccDS;
            // DataGraph *sub = new DataGraph();
            // if (rel->first == 1 && curK == 3)
            // {
            //     bool flag = true;
            // }
            subG->getCC(ccDS, vDS);
            // subG->getSubG(ccDS, DS, maxRV);
            int kMax = 0;
            // if (curK == 15)
            // {
            //     bool flag = true;
            // }
            if (ccDS.size() == 1)
            {
                // cout << "Only one connected component" << endl;
                unordered_map<int, unordered_map<int, int>> DS; // the densest subgraph
                subG->support(DS, vDS);
                // vector<int> vertices(vDS.begin(), vDS.end());
                // DataGraph *sub = subG->getSubG(vertices);
                // TrussDecomposition *td = new TrussDecomposition(sub);
                double maxRho = rhoDS;
                set<int> maxRV = vDS;
                // set<int> maxR = RDS;
                // unordered_map<int, double> maxVert2weight = v2wDS;
                // get kMax-truss holding maxRho or larger rho
                GKMDS(kMax, maxRV, maxRho, RDS, v2wDS, simG, DS); // RDS, DS and v2wDS不需要返回新值，所以可以不单独赋值，传入函数的是副本
                SkyGroupCand c = SkyGroupCand();
                if (kMax < curK + 1)
                {
                    // cout << "curK: " << curK << endl;
                    expand(c, rhoDS, vDS, RDS, simG, DS, subG, RCCSubGDS, V2WCCSubGDS, v2wDS, curK + 1, kMax);
                    // output kMax, maxRho, maxRV
                    c.rel = rel->first;
                    // cand.push_back(c);
                }
                else
                {
                    // is result, output kMax, maxRho, maxRV

                    c.rel = rel->first;
                    c.k = kMax;
                    c.rho = maxRho;
                    for (int v : maxRV)
                    {
                        c.vertices.insert(subG->seq2id[v]);
                    }
                }
                curK = c.k;
                cand.push_back(c);
            }
            else
            {
                // cout << "Some connected components" << endl;
                // 连通分量可能不包含rel点，所以从包含rel的分量中找到rho最大的分量，判断最大k值是否大于k，如果不是则拓展。然后从剩余的分量里找最大，权重最短路径连接。
                // 用2-近似得到密度最大，而不用k-truss密度最大是因为k-truss不是单调的，通过剥蒜式缩减得到的密度可能不是最大的，即其子图中可能还包含一个密度更大的k-truss。
                // unordered_map<int, double> ccDs2Rho;
                vector<pair<int, double>> ccDs2Rho;
                unordered_map<int, unordered_map<int, double>> ccDs2vert2weight;
                bool isSingle = true;
                for (auto c : ccDS)
                {
                    if (c.second.size() > 1)
                        isSingle = false;
                }
                if (isSingle)
                {
                    for (auto c : ccDS)
                    {
                        ccDs2vert2weight.emplace(c.first, unordered_map<int, double>{{*c.second.begin(), 0}});
                        ccDs2Rho.push_back(make_pair(c.first, 0));
                    }
                    sort(ccDs2Rho.begin(), ccDs2Rho.end(), [&v2wDS](const pair<int, double> &a, const pair<int, double> &b)
                         { return v2wDS[a.first] > v2wDS[b.first]; });
                }
                else
                {
                    for (auto c : ccDS)
                    {
                        ccDs2vert2weight.emplace(c.first, unordered_map<int, double>());
                        auto &vert2weight = ccDs2vert2weight[c.first];
                        // whether contain rel vertices
                        set<int> &curCC = c.second;
                        double total = 0;
                        for (int v : curCC)
                        {
                            double w = 0;
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : simG[v])
                                {
                                    if (!curCC.count(u.first))
                                        continue;
                                    w += u.second + denStyle;
                                }
                            vert2weight.emplace(v, w);
                            total += w;
                        }
                        total /= (2 * curCC.size());
                        ccDs2Rho.push_back(make_pair(c.first, total));
                    }
                    sort(ccDs2Rho.begin(), ccDs2Rho.end(), compareByDensity);
                }

                // find the max cc containing rel
                double aRho = 0;
                set<int> RA;
                int ccRDsmaxRho;
                for (auto &c : ccDs2Rho)
                {
                    int ccId = c.first;
                    set<int> &curCC = ccDS[ccId];
                    set<int> com;
                    set_intersection(curCC.begin(), curCC.end(), RDS.begin(), RDS.end(), inserter(com, com.end()));

                    if (!com.empty())
                    {
                        ccRDsmaxRho = ccId;
                        RA = com;

                        aRho = c.second;
                        break;
                    }
                }
                set<int> &a = ccDS[ccRDsmaxRho];
                unordered_map<int, double> &aW = ccDs2vert2weight[ccRDsmaxRho];
                double maxRho = aRho;
                set<int> maxRV = a;
                // set<int> maxR = ccRDs;
                unordered_map<int, unordered_map<int, int>> aSup;
                SkyGroupCand tempc = SkyGroupCand(); // store the current densest k-truss
                if (!isSingle)
                {

                    // check wether is k-truss if not compute the gap between trussness of vertices in the max cc and the expected k
                    subG->support(aSup, a); // only the edges in cc of maxRV have been counted
                    // 判断是否满足k约束，这里传入的是副本，不对原变量进行修改，如果满足，返回k值、子图、密度
                    // vector<int> vertices(a.begin(), a.end());
                    // DataGraph *sub = subG->getSubG(vertices);
                    // TrussDecomposition *td = new TrussDecomposition(sub);
                    // for (auto &kt : td->k2edge)
                    // {
                    //     if (kt.second.size() != 0)
                    //     {
                    //         cout << "k value: " << kt.first << endl;
                    //         break;
                    //     }
                    // }
                    GKMDS(kMax, maxRV, maxRho, RA, aW, simG, aSup);
                    // SkyGroupCand tempc = SkyGroupCand(); // store the current densest k-truss
                    if (kMax < curK + 1)
                    {
                        // double expRho = aRho;
                        // set<int> expRV = a;
                        // set<int> expR = RA;
                        // unordered_map<int, unordered_map<int, int>> expSup = aSup;
                        // unordered_map<int, double> expW = aW;
                        // expand(tempc, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, curK + 1, kMax); // cannot change a
                        expand(tempc, aRho, a, RA, simG, aSup, subG, RCCSubGDS, V2WCCSubGDS, aW, curK + 1, kMax);
                    }
                    else
                    {
                        // is a candidate
                        // c.rel = rel->first;
                        tempc.k = kMax;
                        tempc.rho = maxRho;
                        for (int v : maxRV)
                        {
                            tempc.vertices.insert(subG->seq2id[v]);
                        }
                        // cand.push_back(c);
                    }
                }
                for (auto &c : ccDs2Rho)
                {
                    if (c.first == ccRDsmaxRho)
                        continue;
                    // use least high weight vertices in cc of subG to connect the current cc with the max cc in the rest of ccs, and loop again.
                    // first check whether the expend a connects with b, if so, do not need to connect again.
                    set<int> &b = ccDS[c.first];
                    set<int> bRest;
                    set_difference(b.begin(), b.end(), a.begin(), a.end(), inserter(bRest, bRest.end()));
                    if (bRest.empty())
                        continue;
                    else if (bRest.size() == b.size())
                    {
                        unordered_map<int, unordered_map<int, double>> transG;
                        vector<int> path;
                        int src = -1, dst = -2;

                        bool connt = constructGraph(src, dst, ccIdSubG, subG, V2WCCSubGDS, transG, a, b); // construct graph on the subgraph containing DS
                        if (!connt)
                            dijkstra(src, dst, transG, path); // get added vertices
                        // cout << "path: " << path.size() << endl;
                        // calculate weight of new connected component
                        unordered_map<int, double> &bW = ccDs2vert2weight[c.first];
                        for (auto &vw : bW)
                        {
                            aW.emplace(vw);
                        }
                        aRho = aRho * a.size() + c.second * b.size();
                        for (int v : b)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!a.count(u.first))
                                        continue;
                                    aW[u.first] += u.second + denStyle;
                                    aW[v] += u.second + denStyle;
                                    aRho += u.second + denStyle;
                                }
                        }
                        a.insert(b.begin(), b.end());

                        double tt = 0;
                        // weight sum of shortest path
                        for (int v : path)
                        {
                            double w = 0;
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!count(path.begin(), path.end(), u.first))
                                        continue;
                                    w += u.second + denStyle;
                                }
                            aW.emplace(v, w);
                            tt += w;
                        }
                        aRho += tt / 2;
                        for (int v : path)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!a.count(u.first))
                                        continue;
                                    aW[u.first] += u.second + denStyle;
                                    aW[v] += u.second + denStyle;
                                    aRho += u.second + denStyle;
                                }
                        }
                        for (int v : path)
                        {
                            a.insert(v);
                        }
                        aRho /= a.size();

                        // update R
                        for (int v : RCCSubGDS)
                        {
                            if (count(path.begin(), path.end(), v))
                                RA.insert(v);
                            if (b.count(v))
                                RA.insert(v);
                        }
                        b.clear();
                    }
                    else // already connected
                    {

                        for (int v : bRest)
                        {
                            aW.emplace(v, 0);
                        }
                        aRho = aRho * a.size();
                        for (int v : bRest)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &U : it->second)
                                {
                                    int u = U.first;
                                    if (a.count(u))
                                    {
                                        aW[u] += U.second + denStyle;
                                        aW[v] += U.second + denStyle;
                                        aRho += U.second + denStyle;
                                    }
                                    else if (v < u && bRest.count(u))
                                    {
                                        aW[v] += U.second + denStyle;
                                        aW[u] += U.second + denStyle;
                                        aRho += U.second + denStyle;
                                    }
                                }
                        }
                        a.insert(bRest.begin(), bRest.end());
                        aRho /= a.size();

                        // update R
                        for (int v : RCCSubGDS)
                        {
                            if (b.count(v))
                                RA.insert(v);
                        }
                    }
                    maxRho = aRho;
                    maxRV = a;
                    // maxVert2weight = aW;
                    // update support
                    aSup.clear();
                    subG->support(aSup, a); // only the edges in cc of maxRV have been counted
                    // vector<int> vertices(a.begin(), a.end());
                    // DataGraph *sub = subG->getSubG(vertices);
                    // TrussDecomposition *td = new TrussDecomposition(sub);
                    // for (auto &kt : td->k2edge)
                    // {
                    //     if (kt.second.size() != 0)
                    //     {
                    //         cout << "k value: " << kt.first << endl;
                    //         break;
                    //     }
                    // }
                    GKMDS(kMax, maxRV, maxRho, RA, aW, simG, aSup);
                    // expend the max cc.
                    if (kMax < curK + 1)
                    {
                        // double expRho = aRho;
                        // set<int> expRV = a;
                        // set<int> expR = RA;
                        // unordered_map<int, unordered_map<int, int>> expSup = aSup;
                        // unordered_map<int, double> expW = aW;
                        // SkyGroupCand expCand = SkyGroupCand();
                        // expand(expCand, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, curK + 1, kMax);
                        // if (expCand.rho >= tempc.rho)
                        // {
                        //     tempc.k = expCand.k;
                        //     tempc.rho = expCand.rho;
                        //     tempc.vertices.clear();
                        //     for (int v : expCand.vertices)
                        //     {
                        //         tempc.vertices.insert(subG->seq2id[v]);
                        //     }
                        // }
                        SkyGroupCand expCand = SkyGroupCand();
                        expand(expCand, aRho, a, RA, simG, aSup, subG, RCCSubGDS, V2WCCSubGDS, aW, curK + 1, kMax);
                        if (expCand.rho >= tempc.rho)
                        {
                            tempc.k = expCand.k;
                            tempc.rho = expCand.rho;
                            tempc.vertices = expCand.vertices;
                            // tempc.vertices.clear();
                            // for (int v : expCand.vertices)
                            // {
                            //     tempc.vertices.insert(subG->seq2id[v]);
                            // }
                        }
                    }
                    else
                    {
                        if (maxRho >= tempc.rho)
                        {
                            tempc.k = kMax;
                            tempc.rho = maxRho;
                            tempc.vertices.clear();
                            for (int v : maxRV)
                            {
                                tempc.vertices.insert(subG->seq2id[v]);
                            }
                        }
                    }
                }
                tempc.rel = rel->first;
                cand.push_back(tempc);
                // cout << "k: " << curK << ", K: " << tempc.k << endl;
                curK = tempc.k;
            }
        }

        // cout << "Find candidate skyline groups" << endl;
        
        skyline(result, cand);
        for (const auto &point : cand)
        {
            cout << "rel: " << point.rel << " k: " << point.k << " rho: " << point.rho << " maxRV: ";
            for (auto v : point.vertices)
                cout << v << ", ";
            cout << endl;
        }

        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "skyline2 : " << duration.count() << " milliseconds." << std::endl;
        ++rel;
        delete kt;
    }

    delete relSubG;
}
void Best1(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
    auto rel = rel2edge.begin();
    DataGraph *relSubG = new DataGraph();

    ///////////////////////////////////////
    while (rel != rel2edge.end())
    {
        for (auto &e : rel->second)
        {
            int src = e.first;
            int dst = e.second;
            relSubG->addEdgeNoMatinC(src, dst);
        }

        TrussDecomposition *kt = new TrussDecomposition(relSubG);

        int curK = k;
        // find k+1-truss

        DataGraph *subG = new DataGraph();
        unordered_map<int, unordered_map<int, double>> simG;
        for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
        {
            int k = it->first;
            if (k <= curK)
            {
                it = kt->k2edge.erase(it);
            }
            else
            {
                for (auto &e : it->second)
                {
                    int src = relSubG->seq2id[e.first];
                    int dst = relSubG->seq2id[e.second];
                    subG->addEdgeAndMainConnect(src, dst, simE, simG);
                }

                ++it;
            }
        }
        double rhoDS = 0;
        set<int> vDS;
        set<int> RDS;
        unordered_map<int, double> v2wDS;
        unordered_map<int, double> V2WCCSubGDS;
        set<int> RCCSubGDS;
        int ccIdSubG = 0;
        // find the densest subgraph, get the maxRho, maxRV, maxR, maxVert2weight, gCCId
        GDS(rhoDS, vDS, RDS, v2wDS, V2WCCSubGDS, ccIdSubG, RCCSubGDS, simG, subG, rel, rel2verts);
        if (vDS.empty())
            break;
        // gccid record the id of the connected component from densest subgraph
        // create a new subgraph
        unordered_map<int, set<int>> ccDS;

        subG->getCC(ccDS, vDS);
        // subG->getSubG(ccDS, DS, maxRV);
        int kMax = 0;

        if (ccDS.size() == 1)
        {
            // cout << "Only one connected component" << endl;
            unordered_map<int, unordered_map<int, int>> DS; // the densest subgraph
            subG->support(DS, vDS);

            double maxRho = rhoDS;
            set<int> maxRV = vDS;

            // get kMax-truss holding maxRho or larger rho
            GKMDS(kMax, maxRV, maxRho, RDS, v2wDS, simG, DS); // RDS, DS and v2wDS不需要返回新值，所以可以不单独赋值，传入函数的是副本
            SkyGroupCand c = SkyGroupCand();
            if (kMax < curK + 1)
            {
                // cout << "curK: " << curK << endl;
                expand(c, rhoDS, vDS, RDS, simG, DS, subG, RCCSubGDS, V2WCCSubGDS, v2wDS, curK + 1, kMax);
                // output kMax, maxRho, maxRV
                c.rel = rel->first;
                // cand.push_back(c);
            }
            else
            {
                // is result, output kMax, maxRho, maxRV

                c.rel = rel->first;
                c.k = kMax;
                c.rho = maxRho;
                for (int v : maxRV)
                {
                    c.vertices.insert(subG->seq2id[v]);
                }
            }
            curK = c.k;
            result.push_back(c);
        }
        else
        {
            // cout << "Some connected components" << endl;
            // 连通分量可能不包含rel点，所以从包含rel的分量中找到rho最大的分量，判断最大k值是否大于k，如果不是则拓展。然后从剩余的分量里找最大，权重最短路径连接。
            // 用2-近似得到密度最大，而不用k-truss密度最大是因为k-truss不是单调的，通过剥蒜式缩减得到的密度可能不是最大的，即其子图中可能还包含一个密度更大的k-truss。
            // unordered_map<int, double> ccDs2Rho;
            vector<pair<int, double>> ccDs2Rho;
                unordered_map<int, unordered_map<int, double>> ccDs2vert2weight;
                bool isSingle = true;
                for (auto c : ccDS)
                {
                    if (c.second.size() > 1)
                        isSingle = false;
                }
                if (isSingle)
                {
                    for (auto c : ccDS)
                    {
                        ccDs2vert2weight.emplace(c.first, unordered_map<int, double>{{*c.second.begin(), 0}});
                        ccDs2Rho.push_back(make_pair(c.first, 0));
                    }
                    sort(ccDs2Rho.begin(), ccDs2Rho.end(), [&v2wDS](const pair<int, double> &a, const pair<int, double> &b)
                         { return v2wDS[a.first] > v2wDS[b.first]; });
                }
                else
                {
                    for (auto c : ccDS)
                    {
                        ccDs2vert2weight.emplace(c.first, unordered_map<int, double>());
                        auto &vert2weight = ccDs2vert2weight[c.first];
                        // whether contain rel vertices
                        set<int> &curCC = c.second;
                        double total = 0;
                        for (int v : curCC)
                        {
                            double w = 0;
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : simG[v])
                                {
                                    if (!curCC.count(u.first))
                                        continue;
                                    w += u.second + denStyle;
                                }
                            vert2weight.emplace(v, w);
                            total += w;
                        }
                        total /= (2 * curCC.size());
                        ccDs2Rho.push_back(make_pair(c.first, total));
                    }
                    sort(ccDs2Rho.begin(), ccDs2Rho.end(), compareByDensity);
                }

                // find the max cc containing rel
                double aRho = 0;
                set<int> RA;
                int ccRDsmaxRho;
                for (auto &c : ccDs2Rho)
                {
                    int ccId = c.first;
                    set<int> &curCC = ccDS[ccId];
                    set<int> com;
                    set_intersection(curCC.begin(), curCC.end(), RDS.begin(), RDS.end(), inserter(com, com.end()));

                    if (!com.empty())
                    {
                        ccRDsmaxRho = ccId;
                        RA = com;

                        aRho = c.second;
                        break;
                    }
                }
                set<int> &a = ccDS[ccRDsmaxRho];
                unordered_map<int, double> &aW = ccDs2vert2weight[ccRDsmaxRho];
                double maxRho = aRho;
                set<int> maxRV = a;
                // set<int> maxR = ccRDs;
                unordered_map<int, unordered_map<int, int>> aSup;
                SkyGroupCand tempc = SkyGroupCand(); // store the current densest k-truss
                if (!isSingle)
                {

                    // check wether is k-truss if not compute the gap between trussness of vertices in the max cc and the expected k
                    subG->support(aSup, a); // only the edges in cc of maxRV have been counted
                    // 判断是否满足k约束，这里传入的是副本，不对原变量进行修改，如果满足，返回k值、子图、密度
                    // vector<int> vertices(a.begin(), a.end());
                    // DataGraph *sub = subG->getSubG(vertices);
                    // TrussDecomposition *td = new TrussDecomposition(sub);
                    // for (auto &kt : td->k2edge)
                    // {
                    //     if (kt.second.size() != 0)
                    //     {
                    //         cout << "k value: " << kt.first << endl;
                    //         break;
                    //     }
                    // }
                    GKMDS(kMax, maxRV, maxRho, RA, aW, simG, aSup);
                    // SkyGroupCand tempc = SkyGroupCand(); // store the current densest k-truss
                    if (kMax < curK + 1)
                    {
                        // double expRho = aRho;
                        // set<int> expRV = a;
                        // set<int> expR = RA;
                        // unordered_map<int, unordered_map<int, int>> expSup = aSup;
                        // unordered_map<int, double> expW = aW;
                        // expand(tempc, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, curK + 1, kMax); // cannot change a
                        expand(tempc, aRho, a, RA, simG, aSup, subG, RCCSubGDS, V2WCCSubGDS, aW, curK + 1, kMax);
                    }
                    else
                    {
                        // is a candidate
                        // c.rel = rel->first;
                        tempc.k = kMax;
                        tempc.rho = maxRho;
                        for (int v : maxRV)
                        {
                            tempc.vertices.insert(subG->seq2id[v]);
                        }
                        // cand.push_back(c);
                    }
                }
                for (auto &c : ccDs2Rho)
                {
                    if (c.first == ccRDsmaxRho)
                        continue;
                    // use least high weight vertices in cc of subG to connect the current cc with the max cc in the rest of ccs, and loop again.
                    // first check whether the expend a connects with b, if so, do not need to connect again.
                    set<int> &b = ccDS[c.first];
                    set<int> bRest;
                    set_difference(b.begin(), b.end(), a.begin(), a.end(), inserter(bRest, bRest.end()));
                    if (bRest.empty())
                        continue;
                    else if (bRest.size() == b.size())
                    {
                        unordered_map<int, unordered_map<int, double>> transG;
                        vector<int> path;
                        int src = -1, dst = -2;

                        bool connt = constructGraph(src, dst, ccIdSubG, subG, V2WCCSubGDS, transG, a, b); // construct graph on the subgraph containing DS
                        if (!connt)
                            dijkstra(src, dst, transG, path); // get added vertices
                        // cout << "path: " << path.size() << endl;
                        // calculate weight of new connected component
                        unordered_map<int, double> &bW = ccDs2vert2weight[c.first];
                        for (auto &vw : bW)
                        {
                            aW.emplace(vw);
                        }
                        aRho = aRho * a.size() + c.second * b.size();
                        for (int v : b)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!a.count(u.first))
                                        continue;
                                    aW[u.first] += u.second + denStyle;
                                    aW[v] += u.second + denStyle;
                                    aRho += u.second + denStyle;
                                }
                        }
                        a.insert(b.begin(), b.end());

                        double tt = 0;
                        // weight sum of shortest path
                        for (int v : path)
                        {
                            double w = 0;
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!count(path.begin(), path.end(), u.first))
                                        continue;
                                    w += u.second + denStyle;
                                }
                            aW.emplace(v, w);
                            tt += w;
                        }
                        aRho += tt / 2;
                        for (int v : path)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &u : it->second)
                                {
                                    if (!a.count(u.first))
                                        continue;
                                    aW[u.first] += u.second + denStyle;
                                    aW[v] += u.second + denStyle;
                                    aRho += u.second + denStyle;
                                }
                        }
                        for (int v : path)
                        {
                            a.insert(v);
                        }
                        aRho /= a.size();

                        // update R
                        for (int v : RCCSubGDS)
                        {
                            if (count(path.begin(), path.end(), v))
                                RA.insert(v);
                            if (b.count(v))
                                RA.insert(v);
                        }
                        b.clear();
                    }
                    else // already connected
                    {

                        for (int v : bRest)
                        {
                            aW.emplace(v, 0);
                        }
                        aRho = aRho * a.size();
                        for (int v : bRest)
                        {
                            auto it = simG.find(v);
                            if (it != simG.end())
                                for (auto &U : it->second)
                                {
                                    int u = U.first;
                                    if (a.count(u))
                                    {
                                        aW[u] += U.second + denStyle;
                                        aW[v] += U.second + denStyle;
                                        aRho += U.second + denStyle;
                                    }
                                    else if (v < u && bRest.count(u))
                                    {
                                        aW[v] += U.second + denStyle;
                                        aW[u] += U.second + denStyle;
                                        aRho += U.second + denStyle;
                                    }
                                }
                        }
                        a.insert(bRest.begin(), bRest.end());
                        aRho /= a.size();

                        // update R
                        for (int v : RCCSubGDS)
                        {
                            if (b.count(v))
                                RA.insert(v);
                        }
                    }
                    maxRho = aRho;
                    maxRV = a;
                    // maxVert2weight = aW;
                    // update support
                    aSup.clear();
                    subG->support(aSup, a); // only the edges in cc of maxRV have been counted
                    // vector<int> vertices(a.begin(), a.end());
                    // DataGraph *sub = subG->getSubG(vertices);
                    // TrussDecomposition *td = new TrussDecomposition(sub);
                    // for (auto &kt : td->k2edge)
                    // {
                    //     if (kt.second.size() != 0)
                    //     {
                    //         cout << "k value: " << kt.first << endl;
                    //         break;
                    //     }
                    // }
                    GKMDS(kMax, maxRV, maxRho, RA, aW, simG, aSup);
                    // expend the max cc.
                    if (kMax < curK + 1)
                    {
                        // double expRho = aRho;
                        // set<int> expRV = a;
                        // set<int> expR = RA;
                        // unordered_map<int, unordered_map<int, int>> expSup = aSup;
                        // unordered_map<int, double> expW = aW;
                        // SkyGroupCand expCand = SkyGroupCand();
                        // expand(expCand, expRho, expRV, expR, simG, expSup, subG, RCCSubGDS, V2WCCSubGDS, expW, curK + 1, kMax);
                        // if (expCand.rho >= tempc.rho)
                        // {
                        //     tempc.k = expCand.k;
                        //     tempc.rho = expCand.rho;
                        //     tempc.vertices.clear();
                        //     for (int v : expCand.vertices)
                        //     {
                        //         tempc.vertices.insert(subG->seq2id[v]);
                        //     }
                        // }
                        SkyGroupCand expCand = SkyGroupCand();
                        expand(expCand, aRho, a, RA, simG, aSup, subG, RCCSubGDS, V2WCCSubGDS, aW, curK + 1, kMax);
                        if (expCand.rho >= tempc.rho)
                        {
                            tempc.k = expCand.k;
                            tempc.rho = expCand.rho;
                            tempc.vertices = expCand.vertices;
                            // tempc.vertices.clear();
                            // for (int v : expCand.vertices)
                            // {
                            //     tempc.vertices.insert(subG->seq2id[v]);
                            // }
                        }
                    }
                    else
                    {
                        if (maxRho >= tempc.rho)
                        {
                            tempc.k = kMax;
                            tempc.rho = maxRho;
                            tempc.vertices.clear();
                            for (int v : maxRV)
                            {
                                tempc.vertices.insert(subG->seq2id[v]);
                            }
                        }
                    }
                }
                tempc.rel = rel->first;
                result.push_back(tempc);
                // cout << "k: " << curK << ", K: " << tempc.k << endl;
                curK = tempc.k;
            }
        ++rel;
        delete kt;
    }

    delete relSubG;
}