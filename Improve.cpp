
#include <iostream>
#include <ctime>
#include <string>
#include <set>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include <memory>
#include <limits.h>
// #include "DataGraph.h"
#include "Improve.h"
// void dfs(int node, unordered_set<int> &visited, unordered_map<int, unordered_map<int, int>> &G)
// {
//     visited.insert(node);

//     for (auto &e : G[node])
//     {
//         int neighbor = e.first;
//         if (visited.find(neighbor) == visited.end())
//         {
//             dfs(neighbor, visited, G);
//         }
//     }
// }
// return the maxRho, maxRV, maxR, maxVert2weight
void GDCS(double &maxRho, set<int> &maxRV, set<int> &maxR, unordered_map<int, double> &maxVert2weight, unordered_map<int, unordered_map<int, int>> &maxE2sups,
          unordered_map<int, unordered_map<int, double>> &simG, DataGraph *relSubG, map<int, vector<Edge>, greater<int>>::iterator rel,
          vector<SkyGroupCand> &cand, unordered_map<int, set<int>> &rel2verts, int k, unordered_map<int, unordered_map<int, unordered_map<int, int>>> &e2sups)
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
    for (auto cc : ccId2RelVerts)
    {
        // cout << "Connected component " << i++ << ": " << endl;
        // if(k==3&&i==164){
        //     bool flag = true;
        // }
        unordered_map<int, unordered_map<int, int>> &e2sup = e2sups[cc.first];
        set<int> curCC = relSubG->CC[cc.first];
        set<int> R;
        for (auto &v : cc.second)
            R.insert(relSubG->id2seq[v]);

        /////////////////////////density of whole graph
        int vN = curCC.size();
        unordered_map<int, double> vert2weight;
        set<int> removeV;

        // record weight for each vertices in this connected component;
        double maxWeight = 0;
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
            maxWeight = max(maxWeight, w);
        }
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // cout << "maxWeight:" << maxWeight << ", " << duration.count() << " milliseconds." << endl;
        // unordered_map<int, double> vert2weightCopy = vert2weight;
        if (maxWeight <= 2 * maxRho)
        {
            // cout << maxWeight / 2 << " <= " << maxRho << endl;
            // cout << "first pruning is successful" << endl;
            continue;
        }

        double rho = 0;
        for (auto &v : vert2weight)
            rho += v.second;
        rho /= 2 * vN;
        int z = -1;

        unordered_map<int, vector<int>> bccId2verts;
        set<int> cv;
        start = std::chrono::high_resolution_clock::now();
        unordered_map<int, vector<int>> vert2bccId;
        int bccIndex = 0;
        vector<int> subG(relSubG->CC[cc.first].begin(), relSubG->CC[cc.first].end());
        relSubG->BCC(cv, bccId2verts, subG, removeV, -1);
        bccIndex = bccId2verts.size();

        for (auto &it : bccId2verts)
        {
            for (auto &v : it.second)
            {
                auto tmp = vert2bccId.find(v);
                if (tmp == vert2bccId.end())
                {
                    vert2bccId.emplace(v, vector<int>{it.first});
                }
                else
                    vert2bccId[v].push_back(it.first);
            }
        }

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
        // record the initial rho, vertices and simG, and for find kmax, record the information, containing rho, vertices, rel vertices and weight, of maximum subgraph in the process of deleting vertices
        double tempMaxRho = rho;
        set<int> tempMaxVs = curCC;
        set<int> tempMaxR = R;
        unordered_map<int, unordered_map<int, int>> tmpMaxE2sups = e2sup;
        unordered_map<int, double> tempVert2weight = vert2weight; // in there, delete vertices's weight, but in basic do not delete vertices, just do not update those vertices not in result graph
        // cout << "Rho: " << tempMaxRho << endl;
        bool remove = true;
        // if(i==20){
        //     bool flag = false;
        // }
        // int cnt = 0;
        while (remove) // the termination condition is the pruning condition, so the subgraph size is larger than 2.
        {
            // cnt++;
            // if(k==3&&i==164)cout << "Iteration " << cnt << ": " << endl;
            // if (i == 20 && cnt == 88)
            // {
            //     bool flag = false;
            // }
            // if (k == 15 && cnt == 6)
            // {
            //     bool flag = false;
            // }
            // cout << "cut vertices:" << cv.size() << endl;
            // int vMin = -1;
            remove = false;
            vector<int> canDV;
            set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
            if (canDV.size() == 0)
                break;
            sort(canDV.begin(), canDV.end(), [&](int key1, int key2)
                 { return vert2weight[key1] < vert2weight[key2]; });
            int vWLast = canDV.back();
            maxWeight = vert2weight[vWLast];
            for (int v : cv)
            {
                maxWeight = max(maxWeight, vert2weight[v]);
            }
            if (z != -1)
                maxWeight = max(maxWeight, vert2weight[z]);
            // cout << "maxWeight: " << maxWeight << endl;
            if (maxWeight <= 2 * maxRho)
            {
                // cout << maxWeight / 2 << " <= " << maxRho << endl;
                // cout << "third pruning is successful" << endl;
                break;
            }
            for (int v : canDV)
            {
                // int v = vMin1;
                if (R.size() == 1 && R.count(v))
                {
                    z = v;
                    curCC.erase(v);
                    // v = vMin2;
                }
                else if (checkInd(v, k, e2sup))
                // else if (checkInd(v, relSubG, removeV, k, e2sup))
                {
                    // vMin = v;
                    // if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, v))
                    // {
                    //     cout << "fourth pruning is successful" << endl;
                    //     break;
                    // }
                    remove = true;
                    // if(k == 15) cout<< v << endl;

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
                    // cv.clear();
                    // bccs.clear();
                    start = std::chrono::high_resolution_clock::now();

                    set<int> cv1;
                    unordered_map<int, vector<int>> bccId2verts1;
                    int bccId = vert2bccId[v][0]; // the non-cut vertices only belong to a single bcc

                    vector<int> &bcv = bccId2verts[bccId];
                    sort(bcv.begin(), bcv.end());
                    set<int> bcvSet;
                    set_difference(bcv.begin(), bcv.end(), removeV.begin(), removeV.end(), inserter(bcvSet, bcvSet.begin()));
                    BCC(cv1, bccId2verts1, bcvSet, e2sup);
                    // BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, v,relSubG);
                    // BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV,e2sup);
                    // for (int v : cv1)
                    // {
                    //     cout << v << " ";
                    // }
                    // if (!cv1.empty())
                    //     cout << endl;
                    // cv1.clear();
                    // bccId2verts1.clear();
                    // relSubG->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, v);
                    // for(int v : cv1){
                    //     cout << v << " ";
                    // }
                    // if(!cv1.empty())cout<< endl;
                    updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
                    // cout << "new cut vertices:" << cv1.size() << endl;

                    // kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
                    end = std::chrono::high_resolution_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                    // cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
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
                        tmpMaxE2sups = e2sup;
                    }
                    break;
                }
            }
        }

        // return the maximum from each cc.
        if (tempMaxRho > maxRho)
        {
            // cout << "Found a candidate, rho = " << tempMaxRho << endl;
            maxRho = tempMaxRho;
            maxRV = tempMaxVs;
            maxR = tempMaxR;
            maxVert2weight = tempVert2weight;
            maxE2sups = tmpMaxE2sups;
            // for (auto v : maxRV)
            // 	cout << kTGs->seq2id[v] << ", ";
            // cout << endl;
        }
    }
}
// return kMax, kMaxV, maxRho
void GKMDS(int &kMax, set<int> &kMaxV, set<int> &maxRhoG, double maxRho, set<int> &R,
           unordered_map<int, double> &vert2weight, unordered_map<int, unordered_map<int, double>> &simG, unordered_map<int, unordered_map<int, int>> &edge2sup)
{
    // get the support of each edge in maxRhoG
    // unordered_map<int, unordered_map<int, int>> edge2sup;
    map<int, vector<Edge>> sup2edge;
    for (auto U : edge2sup)
    {
        int u = U.first;
        for (auto V : U.second)
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
        dfs(startV, visited, edge2sup);
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
}

void Improve(map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
    // DataGraph *relSubG;

    auto rel = rel2edge.begin();
    DataGraph *relSubG = new DataGraph();
    // TrussDecomposition *kt = new TrussDecomposition();
    vector<SkyGroupCand> cand;

    ///////////////////////////////////////
    while (rel != rel2edge.end())
    {
        auto start = std::chrono::high_resolution_clock::now();
        // unordered_map<int, set<int>> k2vert;
        cand.clear();
        // vector<Edge> newEs;
        for (auto &e : rel->second)
        {
            int src = e.first;
            int dst = e.second;
            // relSubG->addEdgeAndMainConnect(src, dst, simE, simG);
            relSubG->addEdgeNoMatinC(src, dst); // only to get truss edges
        }

        TrussDecomposition *kt = new TrussDecomposition(relSubG);
        // cout << "TrussMaintance completed" << endl;
        int kMax = 1;
        // int kMax = k;
        while (kMax < kt->kMax)
        {
            DataGraph *subG = new DataGraph();
            unordered_map<int, unordered_map<int, double>> simG;
            for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
            {
                int k = it->first;
                if (k <= kMax)
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
            double maxRho = 0;
            set<int> maxRV;
            set<int> maxR;
            unordered_map<int, double> maxVert2weight;
            unordered_map<int, unordered_map<int, int>> e2sup;
            // get support of each edge
            unordered_map<int, unordered_map<int, unordered_map<int, int>>> e2sups;
            unordered_map<int, vector<Edge>> sup2edge;
            // if (kMax == 14)
            // {
            //     int flag = 0;
            // }
            subG->support(e2sups, sup2edge);
            // if (kMax == 14)
            // {
            //     for (auto &e : e2sups)
            //     {

            //         for (auto &e1 : e.second)
            //             for (auto &e2 : e1.second)
            //                 cout << e1.first << " " << e2.first << " " << e2.second << endl;
            //     }
            // }
            // TrussDecomposition *kt = new TrussDecomposition(subG);
            GDCS(maxRho, maxRV, maxR, maxVert2weight, e2sup, simG, subG, rel, cand, rel2verts, kMax + 1, e2sups);
            if (maxRV.empty())
                break; // dont contain rel vertex
            // GDCS(maxRho, maxRV, maxR, maxVert2weight, e2sup, simG, subG, rel, cand, rel2verts, 2, e2sups);
            // check the result's trussness
            // if (kMax == 14)
            // {
            //     for (auto &e : e2sup)
            //     {

            //         for (auto &e1 : e.second)
            //             cout << e.first << " " << e1.first << " " << e1.second << endl;
            //     }
            // }
            // map<Edge, int> sups;
            // map<int, vector<Edge>> sup2edge1;
            // subG->support(sups, sup2edge1, maxRV);
            // if (kMax == 14)
            // {
            //     for (auto &e : sups)
            //     {
            //         int u = e.first.first;
            //         int v = e.first.second;
            //         int sup = e.second;
            //         if (e2sup[u][v] != sup)
            //             cout << e.first.first << " " << e.first.second << " real " << e.second << " " << e2sup[u][v] << endl;
            //     }
            // }

            // vector<int> ss(maxRV.begin(), maxRV.end());
            // DataGraph *sub = subG->getSubG(ss);
            // // TrussDecomposition *kt1 = new TrussDecomposition(sub);
            // unordered_map<int, unordered_map<int, unordered_map<int, int>>> e2sups1;
            // unordered_map<int, vector<Edge>> sup2edge1;
            // sub->support(e2sups1, sup2edge1);
            set<int> kMaxV;
            GKMDS(kMax, kMaxV, maxRV, maxRho, maxR, maxVert2weight, simG, e2sup);
            SkyGroupCand c = SkyGroupCand();
            c.rel = rel->first;
            c.k = kMax;
            c.rho = maxRho;
            for (int v : kMaxV)
            {
                c.vertices.insert(subG->seq2id[v]);
            }
            // c.vertices = kMaxV;
            cand.push_back(c);
            // if (maxRho <= 0)
            // {
            //     break;
            // }
            // next loop must create a new subgraph, because some edges may be removed
            // k-truss and kMax+1, noted that k may be larger than kMax
            delete subG;
        }
        // cout << "Find candidate skyline groups" << endl;

        // for (const auto &point : cand)
        // {
        //     cout << "rel: " << point.rel << " k: " << point.k << " rho: " << point.rho << " maxRV: ";
        //     for (auto v : point.vertices)
        //         cout << v << ", ";
        //     cout << endl;
        // }
        // cout << "skyline: " << endl;
        skyline(result, cand);
        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "skyline1 : " << duration.count() << " milliseconds." << std::endl;
        ++rel;
        delete kt;
    }

    delete relSubG;
}
void Improve(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
    // DataGraph *relSubG;

    auto rel = rel2edge.begin();
    DataGraph *relSubG = new DataGraph();
    vector<SkyGroupCand> cand;

    ///////////////////////////////////////
    while (rel != rel2edge.end())
    {
        auto start = std::chrono::high_resolution_clock::now();
        // unordered_map<int, set<int>> k2vert;
        cand.clear();
        // vector<Edge> newEs;
        for (auto &e : rel->second)
        {
            int src = e.first;
            int dst = e.second;
            // relSubG->addEdgeAndMainConnect(src, dst, simE, simG);
            relSubG->addEdgeNoMatinC(src, dst); // only to get truss edges
        }

        TrussDecomposition *kt = new TrussDecomposition(relSubG);
        // cout << "TrussMaintance completed" << endl;
        int kMax = k - 1;
        // int kMax = k;
        while (kMax < kt->kMax)
        {
            DataGraph *subG = new DataGraph();
            unordered_map<int, unordered_map<int, double>> simG;
            for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
            {
                int k = it->first;
                if (k <= kMax)
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
            double maxRho = 0;
            set<int> maxRV;
            set<int> maxR;
            unordered_map<int, double> maxVert2weight;
            unordered_map<int, unordered_map<int, int>> e2sup;
            // get support of each edge
            unordered_map<int, unordered_map<int, unordered_map<int, int>>> e2sups;
            unordered_map<int, vector<Edge>> sup2edge;
            // if (kMax == 14)
            // {
            //     int flag = 0;
            // }
            subG->support(e2sups, sup2edge);
            // if (kMax == 14)
            // {
            //     for (auto &e : e2sups)
            //     {

            //         for (auto &e1 : e.second)
            //             for (auto &e2 : e1.second)
            //                 cout << e1.first << " " << e2.first << " " << e2.second << endl;
            //     }
            // }
            // TrussDecomposition *kt = new TrussDecomposition(subG);
            GDCS(maxRho, maxRV, maxR, maxVert2weight, e2sup, simG, subG, rel, cand, rel2verts, kMax + 1, e2sups);
            if (maxRV.empty())
                break; // dont contain rel vertex
            // GDCS(maxRho, maxRV, maxR, maxVert2weight, e2sup, simG, subG, rel, cand, rel2verts, 2, e2sups);
            // check the result's trussness
            // if (kMax == 14)
            // {
            //     for (auto &e : e2sup)
            //     {

            //         for (auto &e1 : e.second)
            //             cout << e.first << " " << e1.first << " " << e1.second << endl;
            //     }
            // }
            // map<Edge, int> sups;
            // map<int, vector<Edge>> sup2edge1;
            // subG->support(sups, sup2edge1, maxRV);
            // if (kMax == 14)
            // {
            //     for (auto &e : sups)
            //     {
            //         int u = e.first.first;
            //         int v = e.first.second;
            //         int sup = e.second;
            //         if (e2sup[u][v] != sup)
            //             cout << e.first.first << " " << e.first.second << " real " << e.second << " " << e2sup[u][v] << endl;
            //     }
            // }

            // vector<int> ss(maxRV.begin(), maxRV.end());
            // DataGraph *sub = subG->getSubG(ss);
            // // TrussDecomposition *kt1 = new TrussDecomposition(sub);
            // unordered_map<int, unordered_map<int, unordered_map<int, int>>> e2sups1;
            // unordered_map<int, vector<Edge>> sup2edge1;
            // sub->support(e2sups1, sup2edge1);
            set<int> kMaxV;
            GKMDS(kMax, kMaxV, maxRV, maxRho, maxR, maxVert2weight, simG, e2sup);
            SkyGroupCand c = SkyGroupCand();
            c.rel = rel->first;
            c.k = kMax;
            c.rho = maxRho;
            for (int v : kMaxV)
            {
                c.vertices.insert(subG->seq2id[v]);
            }
            // c.vertices = kMaxV;
            cand.push_back(c);
            // if (maxRho <= 0)
            // {
            //     break;
            // }
            // next loop must create a new subgraph, because some edges may be removed
            // k-truss and kMax+1, noted that k may be larger than kMax
            delete subG;
        }
        skyline(result, cand);
        for (const auto &point : cand)
        {
            cout << "rel: " << point.rel << " k: " << point.k << " rho: " << point.rho << " maxRV: ";
            for (auto v : point.vertices)
                cout << v << ", ";
            cout << endl;
        }
        ++rel;
        delete kt;
    }

    delete relSubG;
}
void Improve1(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
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
            relSubG->addEdgeNoMatinC(src, dst); // only to get truss edges
        }

        TrussDecomposition *kt = new TrussDecomposition(relSubG);
        int kMax = k;
        DataGraph *subG = new DataGraph();
        unordered_map<int, unordered_map<int, double>> simG;
        for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
        {
            int k = it->first;
            if (k <= kMax)
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
        double maxRho = 0;
        set<int> maxRV;
        set<int> maxR;
        unordered_map<int, double> maxVert2weight;
        unordered_map<int, unordered_map<int, int>> e2sup;
        // get support of each edge
        unordered_map<int, unordered_map<int, unordered_map<int, int>>> e2sups;
        unordered_map<int, vector<Edge>> sup2edge;

        subG->support(e2sups, sup2edge);
       
        GDCS(maxRho, maxRV, maxR, maxVert2weight, e2sup, simG, subG, rel, result, rel2verts, kMax + 1, e2sups);
        if (maxRV.empty())
            break; // dont contain rel vertex
        
        set<int> kMaxV;
        GKMDS(kMax, kMaxV, maxRV, maxRho, maxR, maxVert2weight, simG, e2sup);
        SkyGroupCand c = SkyGroupCand();
        c.rel = rel->first;
        c.k = kMax;
        c.rho = maxRho;
        for (int v : kMaxV)
        {
            c.vertices.insert(subG->seq2id[v]);
        }
        // c.vertices = kMaxV;
        result.push_back(c);
        
        // next loop must create a new subgraph, because some edges may be removed
        // k-truss and kMax+1, noted that k may be larger than kMax
        delete subG;
        ++rel;
        delete kt;
    }

    delete relSubG;
}