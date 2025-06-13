
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

#include "Basic.h"
void relCand(int kMax, TrussDecomposition *kt, DataGraph *relSubG, map<int, vector<Edge>, greater<int>>::iterator rel, vector<SkyGroupCand> &cand, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
	int r = rel->first;
	double maxRho = 0;

	DataGraph *kTGs = new DataGraph();
	unordered_map<int, unordered_map<int, double>> simG;

	for (int k = kt->kMax; k >= 2; k--)
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto it = kt->k2edge.find(k);
		if (it != kt->k2edge.end())
		{
			for (auto &e : kt->k2edge[k])
			{
				int src = relSubG->seq2id[e.first];
				int dst = relSubG->seq2id[e.second];
				kTGs->addEdgeAndMainConnect(src, dst, simE, simG); // must rebuild the graph because some edges do not exist in the induced graph by vertices
			}
		}
		if (k > kMax)
			continue;
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		// cout << k << "Subgraph construction completed. " << duration.count() << " milliseconds." << endl;
		unordered_map<int, set<int>> ccId2RelVerts;

		start = std::chrono::high_resolution_clock::now();
		for (auto cc : kTGs->CC)
		{
			// 与每个连通分量交
			set<int> verts;
			for (int v : cc.second)
			{
				verts.insert(kTGs->seq2id[v]);
			}
			set<int> comVerts;
			set_intersection(rel2verts[r].begin(), rel2verts[r].end(), verts.begin(), verts.end(), inserter(comVerts, comVerts.begin()));
			if (!comVerts.empty())
			{
				ccId2RelVerts.emplace(cc.first, comVerts);
			}
		}
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		// cout << "Found the connected components containing the relation. " << duration.count() << " milliseconds." << endl;
		// connected components
		double maxR = 0;
		set<int> maxRV;
		int i = 0;
		int canNum = 0;
		if (k == 2)
		{
			bool f = true;
		}
		// for a connected component
		for (auto cc : ccId2RelVerts)
		{
			// cout << "Connected component " << i++ << ": " << endl;
			auto start1 = std::chrono::high_resolution_clock::now();
			set<int> curCC = kTGs->CC[cc.first];
			set<int> R;
			for (auto &v : cc.second)
				R.insert(kTGs->id2seq[v]);
			// map<int, unordered_map<int, double>> simG;
			// kTGs->getSimGs(simThreshold, simG, id2att, cc.first); // vertices not id, is seq
			/////////////////////////density of whole graph
			int vN = curCC.size();
			unordered_map<int, double> vert2weight;
			set<int> removeV;

			// queue<int> lessMRhoV;
			double maxWeight = 0;

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

			// cout << "maxWeight:" << maxWeight << ", " << duration.count() << " milliseconds." << endl;
			unordered_map<int, double> vert2weightCopy = vert2weight;
			if (maxWeight <= 2 * maxRho)
			{
				// cout << maxWeight / 2 << " <= " << maxRho << endl;
				// cout << "maxWeight pruning is successful" << endl;
				continue;
			}
			// if (firstMaxRhoLesser(vert2weightCopy, simG, maxRho))
			// {
			// 	// cout << maxWeight / 2 << " <= " << maxRho << endl;
			// 	cout << "first pruning is successful" << endl;
			// 	continue;
			// }
			// unordered_map<int, unordered_map<int, int>> shrinkG;
			// firstMaxRhoLesser(shrinkG, k, cc.first, kTGs, vert2weightCopy, simG, maxRho);

			double rho = 0;
			for (auto &v : vert2weight)
				rho += v.second;
			rho /= 2 * vN;
			int z = -1;
			// find the first k trussness
			if (rho > maxRho)
			{
				maxRho = rho;
				maxRV = curCC;
			}
			map<Edge, int> edge2sup;
			int minSup = INT_MAX;
			vector<Edge> kSupEs;
			kTGs->support(kSupEs, minSup, edge2sup, cc.first, removeV, k); // must recompute because the triangle number of edges in subgraph is lesser than the original graph
			unordered_map<int, vector<int>> bccId2verts;
			set<int> cv;
			// start = std::chrono::high_resolution_clock::now();
			unordered_map<int, vector<int>> vert2bccId;
			int bccIndex = 0;
			// vector<int> sg = vector<int>(kTGs->CC[cc.first].begin(), kTGs->CC[cc.first].end());
			if (minSup + 2 != vN)
			{
				vector<int> subG(kTGs->CC[cc.first].begin(), kTGs->CC[cc.first].end());
				kTGs->BCC(cv, bccId2verts, subG, removeV, -1);
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
			}
			// for (int i : cv)
			// 	cout << i << " ";
			// cout << endl;
			// end = std::chrono::high_resolution_clock::now();
			// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
			double total = 0;
			bool continue_flag = false;
			int itNum = 0;
			while (minSup != k - 2)
			{
				itNum++;
				// cout << "cut vertices:" << cv.size() << endl;
				vector<int> canDV;
				set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
				// auto t = find(canDV.begin(), canDV.end(), z);
				// if (t != canDV.end())
				// 	canDV.erase(t);
				// double vMin = numeric_limits<double>::max();
				// for (int i : canDV)
				// {
				// 	vMin = min(vMin, vert2weight[i]);
				// }
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
					continue_flag = true;
					break;
				}
				auto it = canDV.begin();
				int vMin = *it;
				if (R.size() == 1 && R.count(vMin))
				{
					z = vMin;
					curCC.erase(vMin);
					vMin = *++it;
				}
				if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, vMin))
				{
					continue_flag = true;
					break;
				}
				// if(it==canDV.end())
				// update support
				sort(kTGs->AdjList[vMin].begin(), kTGs->AdjList[vMin].end());
				for (int u : kTGs->AdjList[vMin])
				{
					if (removeV.count(u))
						continue;
					// Edge e = make_pair(min(vMin, u), max(vMin, u));
					// edge2sup.erase(e);
					sort(kTGs->AdjList[u].begin(), kTGs->AdjList[u].end());
					set<int> WT;
					set_intersection(kTGs->AdjList[vMin].begin(), kTGs->AdjList[vMin].end(), kTGs->AdjList[u].begin(), kTGs->AdjList[u].end(), inserter(WT, WT.begin()));
					for (int w : WT)
					{
						if (w < u || removeV.count(w))
							continue;
						Edge e2 = make_pair(u, w);
						// if (removeV.count(w))
						// 	continue;

						// Edge e2 = make_pair(min(u, w), max(u, w));
						auto it = edge2sup.find(e2);
						if (it != edge2sup.end())
						{
							if (--it->second < minSup)
							{
								minSup = it->second;
							}
							if (it->second == k - 2)
								kSupEs.push_back(e2);
						}
					}
				}
				removeV.insert(vMin);
				curCC.erase(vMin);
				// update weight
				// total += vert2weight[vMin];
				rho = (rho * vN - vert2weight[vMin]) / --vN;
				for (auto &e : simG[vMin])
				{
					if (!curCC.count(e.first) && z != e.first)
						continue;
					// simG[e.first].erase(vMin);
					double t = vert2weight[e.first];

					vert2weight[e.first] = t - (e.second + denStyle);
				}

				// vert2weight.erase(vMin);
				// simG.erase(vMin);

				// cv.clear();
				// bccs.clear();
				// start = std::chrono::high_resolution_clock::now();

				if (minSup + 2 != vN)
				{
					set<int> cv1;
					unordered_map<int, vector<int>> bccId2verts1;
					int bccId = vert2bccId[vMin][0]; // the non-cut vertices only belong to a single bcc
					// kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
					kTGs->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, vMin);
					updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
					// cout << "new cut vertices:" << cv1.size() << endl;
					// for (int i : cv)
					// 	cout << i << " ";
					// cout << endl;
				}

				// end = std::chrono::high_resolution_clock::now();
				// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
				R.erase(vMin);
			}
			if (continue_flag)
			{
				// cout << "second pruning is successful" << endl;
				continue;
			}

			// the first density of k trussness
			// rho = (rho * vN - total) / (vN - removeV.size());
			// vN -= removeV.size();
			// maxRho = max(maxRho, rho);
			if (rho > maxRho)
			{
				maxRho = rho;
				// maxRho = max(maxRho, rho);
				maxRV = curCC;
				if (z != -1 && !curCC.count(z))
					maxRV.insert(z);
			}
			// double tempMaxRho = rho;
			// set<int> tempMaxVs = curCC;
			// if (z != -1)
			// 	tempMaxVs.insert(z);
			// check whether the density is enough
			bool remove = true;
			while (remove)
			{
				itNum++;
				// cout << "cut vertices:" << cv.size() << endl;
				// int vMin = -1;
				remove = false;
				vector<int> canDV;
				set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
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
					if (R.size() == 1 && R.count(v))
					{
						z = v;
						curCC.erase(v);
					}
					else if (checkInd(v, kTGs, removeV, k, edge2sup, kSupEs))
					{
						// vMin = v;
						if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, v))
						{
							// cout << "fourth pruning is successful" << endl;
							break;
						}
						remove = true;

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
						// cv.clear();
						// bccs.clear();
						// start = std::chrono::high_resolution_clock::now();
						if (vN != k)
						{
							set<int> cv1;
							unordered_map<int, vector<int>> bccId2verts1;
							int bccId = vert2bccId[v][0]; // the non-cut vertices only belong to a single bcc
							kTGs->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, v);
							updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
							// cout << "new cut vertices:" << cv1.size() << endl;
							// for (int i : cv)
							// 	cout << i << " ";
							// cout << endl;
						}

						// kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
						// end = std::chrono::high_resolution_clock::now();
						// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
						// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
						R.erase(v);

						// if (rho > tempMaxRho)
						// {
						// 	tempMaxRho = rho;
						// 	// maxRho = max(maxRho, rho);
						// 	tempMaxVs = curCC;
						// 	if (z != -1 && !curCC.count(z))
						// 		tempMaxVs.insert(z);
						// }
						if (rho > maxRho)
						{
							maxRho = rho;
							// maxRho = max(maxRho, rho);
							maxRV = curCC;
							if (z != -1 && !curCC.count(z))
								maxRV.insert(z);
						}

						break;
					}
				}
			}

			// cout << ++canNum << " ";
			// if (tempMaxRho > maxRho)
			// {
			// 	cout << "Found a candidate, rho = " << tempMaxRho << endl;
			// 	maxRho = tempMaxRho;
			// 	maxRV = tempMaxVs;
			// 	// for (auto v : maxRV)
			// 	// 	cout << kTGs->seq2id[v] << ", ";
			// 	// cout << endl;
			// }
			auto end1 = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
			// cout << "time: " << duration.count() << " iterate Num: " << itNum << endl;
		}
		// unordered_map<int, map<int, unordered_map<int, double>>> simGs = kTGs->getSimGs(id2att);

		if (maxRV.size())
		{
			// maxRho = maxR;
			SkyGroupCand c = SkyGroupCand();
			c.rel = r;
			c.k = k;
			c.rho = maxRho;
			c.vertices = maxRV;
			cand.push_back(c);
		}
	}
	delete kTGs;
}
void relCand(int mink, int kMax, TrussDecomposition *kt, DataGraph *relSubG, map<int, vector<Edge>, greater<int>>::iterator rel, vector<SkyGroupCand> &cand, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
	int r = rel->first;
	double maxRho = 0;

	DataGraph *kTGs = new DataGraph();
	unordered_map<int, unordered_map<int, double>> simG;

	for (int k = kt->kMax; k >= mink; k--)
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto it = kt->k2edge.find(k);
		if (it != kt->k2edge.end())
		{
			for (auto &e : kt->k2edge[k])
			{
				int src = relSubG->seq2id[e.first];
				int dst = relSubG->seq2id[e.second];
				kTGs->addEdgeAndMainConnect(src, dst, simE, simG); // must rebuild the graph because some edges do not exist in the induced graph by vertices
			}
		}
		if (k > kMax)
			continue;
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		// cout << k << "Subgraph construction completed. " << duration.count() << " milliseconds." << endl;
		unordered_map<int, set<int>> ccId2RelVerts;

		start = std::chrono::high_resolution_clock::now();
		for (auto cc : kTGs->CC)
		{
			// 与每个连通分量交
			set<int> verts;
			for (int v : cc.second)
			{
				verts.insert(kTGs->seq2id[v]);
			}
			set<int> comVerts;
			set_intersection(rel2verts[r].begin(), rel2verts[r].end(), verts.begin(), verts.end(), inserter(comVerts, comVerts.begin()));
			if (!comVerts.empty())
			{
				ccId2RelVerts.emplace(cc.first, comVerts);
			}
		}
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		// cout << "Found the connected components containing the relation. " << duration.count() << " milliseconds." << endl;
		// connected components
		double maxR = 0;
		set<int> maxRV;
		int i = 0;
		int canNum = 0;
		if (k == 2)
		{
			bool f = true;
		}
		// for a connected component
		for (auto cc : ccId2RelVerts)
		{
			// cout << "Connected component " << i++ << ": " << endl;
			auto start1 = std::chrono::high_resolution_clock::now();
			set<int> curCC = kTGs->CC[cc.first];
			set<int> R;
			for (auto &v : cc.second)
				R.insert(kTGs->id2seq[v]);
			// map<int, unordered_map<int, double>> simG;
			// kTGs->getSimGs(simThreshold, simG, id2att, cc.first); // vertices not id, is seq
			/////////////////////////density of whole graph
			int vN = curCC.size();
			unordered_map<int, double> vert2weight;
			set<int> removeV;

			// queue<int> lessMRhoV;
			double maxWeight = 0;

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

			// cout << "maxWeight:" << maxWeight << ", " << duration.count() << " milliseconds." << endl;
			// unordered_map<int, double> vert2weightCopy = vert2weight;
			if (maxWeight <= 2 * maxRho)
			{
				// cout << maxWeight / 2 << " <= " << maxRho << endl;
				// cout << "maxWeight pruning is successful" << endl;
				continue;
			}
			// if (firstMaxRhoLesser(vert2weightCopy, simG, maxRho))
			// {
			// 	// cout << maxWeight / 2 << " <= " << maxRho << endl;
			// 	cout << "first pruning is successful" << endl;
			// 	continue;
			// }
			// unordered_map<int, unordered_map<int, int>> shrinkG;
			// firstMaxRhoLesser(shrinkG, k, cc.first, kTGs, vert2weightCopy, simG, maxRho);

			double rho = 0;
			for (auto &v : vert2weight)
				rho += v.second;
			rho /= 2 * vN;
			int z = -1;
			// find the first k trussness
			if (rho > maxRho)
			{
				maxRho = rho;
				maxRV = curCC;
			}
			map<Edge, int> edge2sup;
			int minSup = INT_MAX;
			vector<Edge> kSupEs;
			kTGs->support(kSupEs, minSup, edge2sup, cc.first, removeV, k); // must recompute because the triangle number of edges in subgraph is lesser than the original graph
			unordered_map<int, vector<int>> bccId2verts;
			set<int> cv;
			// start = std::chrono::high_resolution_clock::now();
			unordered_map<int, vector<int>> vert2bccId;
			int bccIndex = 0;
			// vector<int> sg = vector<int>(kTGs->CC[cc.first].begin(), kTGs->CC[cc.first].end());
			if (minSup + 2 != vN)
			{
				vector<int> subG(kTGs->CC[cc.first].begin(), kTGs->CC[cc.first].end());
				kTGs->BCC(cv, bccId2verts, subG, removeV, -1);
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
			}
			// for (int i : cv)
			// 	cout << i << " ";
			// cout << endl;
			// end = std::chrono::high_resolution_clock::now();
			// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
			double total = 0;
			bool continue_flag = false;
			int itNum = 0;
			while (minSup != k - 2)
			{
				itNum++;
				// cout << "cut vertices:" << cv.size() << endl;
				vector<int> canDV;
				set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
				// auto t = find(canDV.begin(), canDV.end(), z);
				// if (t != canDV.end())
				// 	canDV.erase(t);
				// double vMin = numeric_limits<double>::max();
				// for (int i : canDV)
				// {
				// 	vMin = min(vMin, vert2weight[i]);
				// }
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
					continue_flag = true;
					break;
				}
				auto it = canDV.begin();
				int vMin = *it;
				if (R.size() == 1 && R.count(vMin))
				{
					z = vMin;
					curCC.erase(vMin);
					vMin = *++it;
				}
				// if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, vMin))
				// {
				// 	continue_flag = true;
				// 	break;
				// }
				// if(it==canDV.end())
				// update support
				sort(kTGs->AdjList[vMin].begin(), kTGs->AdjList[vMin].end());
				for (int u : kTGs->AdjList[vMin])
				{
					if (removeV.count(u))
						continue;
					// Edge e = make_pair(min(vMin, u), max(vMin, u));
					// edge2sup.erase(e);
					sort(kTGs->AdjList[u].begin(), kTGs->AdjList[u].end());
					set<int> WT;
					set_intersection(kTGs->AdjList[vMin].begin(), kTGs->AdjList[vMin].end(), kTGs->AdjList[u].begin(), kTGs->AdjList[u].end(), inserter(WT, WT.begin()));
					for (int w : WT)
					{
						if (w < u || removeV.count(w))
							continue;
						Edge e2 = make_pair(u, w);
						// if (removeV.count(w))
						// 	continue;

						// Edge e2 = make_pair(min(u, w), max(u, w));
						auto it = edge2sup.find(e2);
						if (it != edge2sup.end())
						{
							if (--it->second < minSup)
							{
								minSup = it->second;
							}
							if (it->second == k - 2)
								kSupEs.push_back(e2);
						}
					}
				}
				removeV.insert(vMin);
				curCC.erase(vMin);
				// update weight
				// total += vert2weight[vMin];
				rho = (rho * vN - vert2weight[vMin]) / --vN;
				for (auto &e : simG[vMin])
				{
					if (!curCC.count(e.first) && z != e.first)
						continue;
					// simG[e.first].erase(vMin);
					double t = vert2weight[e.first];

					vert2weight[e.first] = t - (e.second + denStyle);
				}

				// vert2weight.erase(vMin);
				// simG.erase(vMin);

				// cv.clear();
				// bccs.clear();
				// start = std::chrono::high_resolution_clock::now();

				if (minSup + 2 != vN)
				{
					set<int> cv1;
					unordered_map<int, vector<int>> bccId2verts1;
					int bccId = vert2bccId[vMin][0]; // the non-cut vertices only belong to a single bcc
					// kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
					kTGs->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, vMin);
					updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
					// cout << "new cut vertices:" << cv1.size() << endl;
					// for (int i : cv)
					// 	cout << i << " ";
					// cout << endl;
				}

				// end = std::chrono::high_resolution_clock::now();
				// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
				R.erase(vMin);
			}
			// if (continue_flag)
			// {
			// 	// cout << "second pruning is successful" << endl;
			// 	continue;
			// }

			// the first density of k trussness
			// rho = (rho * vN - total) / (vN - removeV.size());
			// vN -= removeV.size();
			// maxRho = max(maxRho, rho);
			if (rho > maxRho)
			{
				maxRho = rho;
				// maxRho = max(maxRho, rho);
				maxRV = curCC;
				if (z != -1 && !curCC.count(z))
					maxRV.insert(z);
			}
			// double tempMaxRho = rho;
			// set<int> tempMaxVs = curCC;
			// if (z != -1)
			// 	tempMaxVs.insert(z);
			// check whether the density is enough
			bool remove = true;
			while (remove)
			{
				itNum++;
				// cout << "cut vertices:" << cv.size() << endl;
				// int vMin = -1;
				remove = false;
				vector<int> canDV;
				set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
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
					if (R.size() == 1 && R.count(v))
					{
						z = v;
						curCC.erase(v);
					}
					else if (checkInd(v, kTGs, removeV, k, edge2sup, kSupEs))
					{
						// vMin = v;
						// if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, v))
						// {
						// 	// cout << "fourth pruning is successful" << endl;
						// 	break;
						// }
						remove = true;

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
						// cv.clear();
						// bccs.clear();
						// start = std::chrono::high_resolution_clock::now();
						if (vN != k)
						{
							set<int> cv1;
							unordered_map<int, vector<int>> bccId2verts1;
							int bccId = vert2bccId[v][0]; // the non-cut vertices only belong to a single bcc
							kTGs->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, v);
							updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
							// cout << "new cut vertices:" << cv1.size() << endl;
							// for (int i : cv)
							// 	cout << i << " ";
							// cout << endl;
						}

						// kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
						// end = std::chrono::high_resolution_clock::now();
						// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
						// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
						R.erase(v);

						// if (rho > tempMaxRho)
						// {
						// 	tempMaxRho = rho;
						// 	// maxRho = max(maxRho, rho);
						// 	tempMaxVs = curCC;
						// 	if (z != -1 && !curCC.count(z))
						// 		tempMaxVs.insert(z);
						// }
						if (rho > maxRho)
						{
							maxRho = rho;
							// maxRho = max(maxRho, rho);
							maxRV = curCC;
							if (z != -1 && !curCC.count(z))
								maxRV.insert(z);
						}

						break;
					}
				}
			}

			// cout << ++canNum << " ";
			// if (tempMaxRho > maxRho)
			// {
			// 	cout << "Found a candidate, rho = " << tempMaxRho << endl;
			// 	maxRho = tempMaxRho;
			// 	maxRV = tempMaxVs;
			// 	// for (auto v : maxRV)
			// 	// 	cout << kTGs->seq2id[v] << ", ";
			// 	// cout << endl;
			// }
			auto end1 = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
			// cout << "time: " << duration.count() << " iterate Num: " << itNum << endl;
		}
		// unordered_map<int, map<int, unordered_map<int, double>>> simGs = kTGs->getSimGs(id2att);

		if (maxRV.size())
		{
			// maxRho = maxR;
			SkyGroupCand c = SkyGroupCand();
			c.rel = r;
			c.k = k;
			c.rho = maxRho;
			c.vertices = maxRV;
			cand.push_back(c);
		}
	}
	delete kTGs;
}
void Basic(map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
	// DataGraph *relSubG;

	auto rel = rel2edge.begin();
	DataGraph *relSubG = new DataGraph();
	// TrussDecomposition *kt = new TrussDecomposition();
	vector<SkyGroupCand> cand;
	///////////////////////////////////////
	while (rel != rel2edge.end())
	{
		// unordered_map<int, set<int>> k2vert;
		cand.clear();
		vector<Edge> newEs;
		relSubG->addEdges(newEs, rel->second);
		// for (auto &e : rel->second)
		// {
		// 	relSubG->addEdge(e);
		// 	kt->UpdateTrussofInsert2(newEs, e, relSubG);
		// }
		TrussDecomposition *kt = new TrussDecomposition(relSubG);
		// cout << "TrussMaintance completed" << endl;

		int kMax = 0; // max trussness of current rel edges
		for (Edge &e : newEs)
		{
			int k = kt->trussd[e];
			kMax = max(kMax, k);
		}
		// cout<<"The maximum trussness of new edges";

		// auto start1 = std::chrono::high_resolution_clock::now();
		relCand(kMax, kt, relSubG, rel, cand, simE, rel2verts);
		// auto end1 = std::chrono::high_resolution_clock::now();
		// auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
		// std::cout << "Find candidates for a relevant : " << duration1.count() << " milliseconds." << std::endl;
		// cout << "Find candidate skyline groups" << endl;
		delete kt;
		skyline(result, cand);
		++rel;
	}
	delete relSubG;
}
void Basic(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
	// DataGraph *relSubG;

	auto rel = rel2edge.begin();
	DataGraph *relSubG = new DataGraph();
	// TrussDecomposition *kt = new TrussDecomposition();
	vector<SkyGroupCand> cand;
	///////////////////////////////////////
	while (rel != rel2edge.end())
	{
		// unordered_map<int, set<int>> k2vert;
		cand.clear();
		vector<Edge> newEs;
		relSubG->addEdges(newEs, rel->second);
		// for (auto &e : rel->second)
		// {
		// 	relSubG->addEdge(e);
		// 	kt->UpdateTrussofInsert2(newEs, e, relSubG);
		// }
		TrussDecomposition *kt = new TrussDecomposition(relSubG);

		// cout << "TrussMaintance completed" << endl;

		int kMax = 0; // max trussness of current rel edges
		for (Edge &e : newEs)
		{
			int tk = kt->trussd[e];
			kMax = max(kMax, tk);
		}
		if (kMax < k)
		{
			++rel;
			continue;
		}

		// cout<<"The maximum trussness of new edges";
		for (auto it = kt->k2edge.begin(); it != kt->k2edge.end();)
		{
			if (it->first < k)
			{
				it = kt->k2edge.erase(it);
			}
			else
			{
				++it;
			}
		}
		// auto start1 = std::chrono::high_resolution_clock::now();
		relCand(k, kMax, kt, relSubG, rel, cand, simE, rel2verts);
		// auto end1 = std::chrono::high_resolution_clock::now();
		// auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
		// std::cout << "Find candidates for a relevant : " << duration1.count() << " milliseconds." << std::endl;
		// cout << "Find candidate skyline groups" << endl;
		delete kt;
		skyline(result, cand);
		for (const auto &point : cand)
        {
            cout << "rel: " << point.rel << " k: " << point.k << " rho: " << point.rho << " maxRV: ";
            for (auto v : point.vertices)
                cout << v << ", ";
            cout << endl;
        }
		++rel;
	}
	delete relSubG;
}
void relCand1(int curk, TrussDecomposition *kt, DataGraph *relSubG, map<int, vector<Edge>, greater<int>>::iterator rel, vector<SkyGroupCand> &cand, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{
	int r = rel->first;
	double maxRho = 0;

	DataGraph *kTGs = new DataGraph();
	unordered_map<int, unordered_map<int, double>> simG;

	for (int k = kt->kMax; k >= curk + 1; k--)
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto it = kt->k2edge.find(k);
		if (it != kt->k2edge.end())
		{
			for (auto &e : kt->k2edge[k])
			{
				int src = relSubG->seq2id[e.first];
				int dst = relSubG->seq2id[e.second];
				kTGs->addEdgeAndMainConnect(src, dst, simE, simG); // must rebuild the graph because some edges do not exist in the induced graph by vertices
			}
		}
		if (k > curk + 1)
			continue;
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		cout << k << "Subgraph construction completed. " << duration.count() << " milliseconds." << endl;
		unordered_map<int, set<int>> ccId2RelVerts;

		start = std::chrono::high_resolution_clock::now();
		for (auto cc : kTGs->CC)
		{
			// 与每个连通分量交
			set<int> verts;
			for (int v : cc.second)
			{
				verts.insert(kTGs->seq2id[v]);
			}
			set<int> comVerts;
			set_intersection(rel2verts[r].begin(), rel2verts[r].end(), verts.begin(), verts.end(), inserter(comVerts, comVerts.begin()));
			if (!comVerts.empty())
			{
				ccId2RelVerts.emplace(cc.first, comVerts);
			}
		}
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		cout << "Found the connected components containing the relation. " << duration.count() << " milliseconds." << endl;
		// connected components
		double maxR = 0;
		set<int> maxRV;
		int i = 0;
		int canNum = 0;

		// for a connected component
		for (auto cc : ccId2RelVerts)
		{
			// cout << "Connected component " << i++ << ": " << endl;
			auto start1 = std::chrono::high_resolution_clock::now();
			set<int> curCC = kTGs->CC[cc.first];
			set<int> R;
			for (auto &v : cc.second)
				R.insert(kTGs->id2seq[v]);
			// map<int, unordered_map<int, double>> simG;
			// kTGs->getSimGs(simThreshold, simG, id2att, cc.first); // vertices not id, is seq
			/////////////////////////density of whole graph
			int vN = curCC.size();
			unordered_map<int, double> vert2weight;
			set<int> removeV;

			// queue<int> lessMRhoV;
			double maxWeight = 0;

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

			// cout << "maxWeight:" << maxWeight << ", " << duration.count() << " milliseconds." << endl;
			unordered_map<int, double> vert2weightCopy = vert2weight;
			if (maxWeight <= 2 * maxRho)
			{
				// cout << maxWeight / 2 << " <= " << maxRho << endl;
				// cout << "maxWeight pruning is successful" << endl;
				continue;
			}
			// if (firstMaxRhoLesser(vert2weightCopy, simG, maxRho))
			// {
			// 	// cout << maxWeight / 2 << " <= " << maxRho << endl;
			// 	cout << "first pruning is successful" << endl;
			// 	continue;
			// }
			// if(i==8){
			// 	bool stop = true;
			// }
			unordered_map<int, unordered_map<int, int>> shrinkG;
			firstMaxRhoLesser(shrinkG, k, cc.first, kTGs, vert2weightCopy, simG, maxRho);

			double rho = 0;
			for (auto &v : vert2weight)
				rho += v.second;
			rho /= 2 * vN;
			int z = -1;
			// find the first k trussness
			if (rho > maxRho)
			{
				maxRho = rho;
				maxRV = curCC;
			}
			map<Edge, int> edge2sup;
			int minSup = INT_MAX;
			vector<Edge> kSupEs;
			kTGs->support(kSupEs, minSup, edge2sup, cc.first, removeV, k); // must recompute because the triangle number of edges in subgraph is lesser than the original graph
			unordered_map<int, vector<int>> bccId2verts;
			set<int> cv;
			// start = std::chrono::high_resolution_clock::now();
			unordered_map<int, vector<int>> vert2bccId;
			int bccIndex = 0;
			// vector<int> sg = vector<int>(kTGs->CC[cc.first].begin(), kTGs->CC[cc.first].end());
			if (minSup + 2 != vN)
			{
				vector<int> subG(kTGs->CC[cc.first].begin(), kTGs->CC[cc.first].end());
				kTGs->BCC(cv, bccId2verts, subG, removeV, -1);
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
			}
			// for (int i : cv)
			// 	cout << i << " ";
			// cout << endl;
			// end = std::chrono::high_resolution_clock::now();
			// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
			double total = 0;
			bool continue_flag = false;
			int itNum = 0;
			while (minSup != k - 2)
			{
				itNum++;
				// cout << "cut vertices:" << cv.size() << endl;
				vector<int> canDV;
				set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
				// auto t = find(canDV.begin(), canDV.end(), z);
				// if (t != canDV.end())
				// 	canDV.erase(t);
				// double vMin = numeric_limits<double>::max();
				// for (int i : canDV)
				// {
				// 	vMin = min(vMin, vert2weight[i]);
				// }
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
					continue_flag = true;
					break;
				}
				auto it = canDV.begin();
				int vMin = *it;
				if (R.size() == 1 && R.count(vMin))
				{
					z = vMin;
					curCC.erase(vMin);
					vMin = *++it;
				}
				if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, vMin))
				{
					continue_flag = true;
					break;
				}
				// if(it==canDV.end())
				// update support
				sort(kTGs->AdjList[vMin].begin(), kTGs->AdjList[vMin].end());
				for (int u : kTGs->AdjList[vMin])
				{
					if (removeV.count(u))
						continue;
					// Edge e = make_pair(min(vMin, u), max(vMin, u));
					// edge2sup.erase(e);
					sort(kTGs->AdjList[u].begin(), kTGs->AdjList[u].end());
					set<int> WT;
					set_intersection(kTGs->AdjList[vMin].begin(), kTGs->AdjList[vMin].end(), kTGs->AdjList[u].begin(), kTGs->AdjList[u].end(), inserter(WT, WT.begin()));
					for (int w : WT)
					{
						if (w < u || removeV.count(w))
							continue;
						Edge e2 = make_pair(u, w);
						// if (removeV.count(w))
						// 	continue;

						// Edge e2 = make_pair(min(u, w), max(u, w));
						auto it = edge2sup.find(e2);
						if (it != edge2sup.end())
						{
							if (--it->second < minSup)
							{
								minSup = it->second;
							}
							if (it->second == k - 2)
								kSupEs.push_back(e2);
						}
					}
				}
				removeV.insert(vMin);
				curCC.erase(vMin);
				// update weight
				// total += vert2weight[vMin];
				rho = (rho * vN - vert2weight[vMin]) / --vN;
				for (auto &e : simG[vMin])
				{
					if (!curCC.count(e.first) && z != e.first)
						continue;
					// simG[e.first].erase(vMin);
					double t = vert2weight[e.first];

					vert2weight[e.first] = t - (e.second + denStyle);
				}

				// vert2weight.erase(vMin);
				// simG.erase(vMin);

				// cv.clear();
				// bccs.clear();
				// start = std::chrono::high_resolution_clock::now();

				if (minSup + 2 != vN)
				{
					set<int> cv1;
					unordered_map<int, vector<int>> bccId2verts1;
					int bccId = vert2bccId[vMin][0]; // the non-cut vertices only belong to a single bcc
					// kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
					kTGs->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, vMin);
					updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
					// cout << "new cut vertices:" << cv1.size() << endl;
					// for (int i : cv)
					// 	cout << i << " ";
					// cout << endl;
				}

				// end = std::chrono::high_resolution_clock::now();
				// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
				R.erase(vMin);
			}
			if (continue_flag)
			{
				// cout << "second pruning is successful" << endl;
				continue;
			}

			// the first density of k trussness
			// rho = (rho * vN - total) / (vN - removeV.size());
			// vN -= removeV.size();
			// maxRho = max(maxRho, rho);
			if (rho > maxRho)
			{
				maxRho = rho;
				// maxRho = max(maxRho, rho);
				maxRV = curCC;
				if (z != -1 && !curCC.count(z))
					maxRV.insert(z);
			}
			// double tempMaxRho = rho;
			// set<int> tempMaxVs = curCC;
			// if (z != -1)
			// 	tempMaxVs.insert(z);
			// check whether the density is enough
			bool remove = true;
			while (remove)
			{
				itNum++;
				// cout << "cut vertices:" << cv.size() << endl;
				// int vMin = -1;
				remove = false;
				vector<int> canDV;
				set_difference(curCC.begin(), curCC.end(), cv.begin(), cv.end(), back_inserter(canDV));
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
					if (R.size() == 1 && R.count(v))
					{
						z = v;
						curCC.erase(v);
					}
					else if (checkInd(v, kTGs, removeV, k, edge2sup, kSupEs))
					{
						// vMin = v;
						if (secondMaxRhoLesser(vert2weightCopy, simG, maxRho, v))
						{
							// cout << "fourth pruning is successful" << endl;
							break;
						}
						remove = true;

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
						// cv.clear();
						// bccs.clear();
						// start = std::chrono::high_resolution_clock::now();
						if (vN != k)
						{
							set<int> cv1;
							unordered_map<int, vector<int>> bccId2verts1;
							int bccId = vert2bccId[v][0]; // the non-cut vertices only belong to a single bcc
							kTGs->BCC(cv1, bccId2verts1, bccId2verts[bccId], removeV, v);
							updateBcc(cv, bccId2verts, cv1, bccId2verts1, vert2bccId, bccId, bccIndex);
							// cout << "new cut vertices:" << cv1.size() << endl;
							// for (int i : cv)
							// 	cout << i << " ";
							// cout << endl;
						}

						// kTGs->BCC(cv, bccs, cc.first, removeV); // have repeated computation
						// end = std::chrono::high_resolution_clock::now();
						// duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
						// cout << "BCC search time: " << duration.count() << " milliseconds." << endl;
						R.erase(v);

						// if (rho > tempMaxRho)
						// {
						// 	tempMaxRho = rho;
						// 	// maxRho = max(maxRho, rho);
						// 	tempMaxVs = curCC;
						// 	if (z != -1 && !curCC.count(z))
						// 		tempMaxVs.insert(z);
						// }
						if (rho > maxRho)
						{
							maxRho = rho;
							// maxRho = max(maxRho, rho);
							maxRV = curCC;
							if (z != -1 && !curCC.count(z))
								maxRV.insert(z);
						}

						break;
					}
				}
			}

			// cout << ++canNum << " ";
			// if (tempMaxRho > maxRho)
			// {
			// 	cout << "Found a candidate, rho = " << tempMaxRho << endl;
			// 	maxRho = tempMaxRho;
			// 	maxRV = tempMaxVs;
			// 	// for (auto v : maxRV)
			// 	// 	cout << kTGs->seq2id[v] << ", ";
			// 	// cout << endl;
			// }
			auto end1 = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
			// cout << "time: " << duration.count() << " iterate Num: " << itNum << endl;
		}
		// unordered_map<int, map<int, unordered_map<int, double>>> simGs = kTGs->getSimGs(id2att);

		if (maxRV.size())
		{
			// maxRho = maxR;
			SkyGroupCand c = SkyGroupCand();
			c.rel = r;
			c.k = k;
			c.rho = maxRho;
			c.vertices = maxRV;
			cand.push_back(c);
		}
	}
	delete kTGs;
}
void Basic1(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts)
{

	auto rel = rel2edge.begin();
	DataGraph *relSubG = new DataGraph();
	///////////////////////////////////////
	while (rel != rel2edge.end())
	{
		vector<Edge> newEs;
		relSubG->addEdges(newEs, rel->second);
		TrussDecomposition *kt = new TrussDecomposition(relSubG);
		relCand1(k, kt, relSubG, rel, result, simE, rel2verts);

		delete kt;
		++rel;
	}
	delete relSubG;
}