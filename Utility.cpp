#include "Utility.h"
#include <queue>
// bool checkInd(int v, DataGraph *g, set<int> &removeV, int k, map<Edge, int> &edge2sup)
// {
// 	vector<Edge> DEdges;
// 	vector<Edge> newSupEs;
// 	for (int u : g->AdjList[v])
// 	{
// 		if (removeV.count(u))
// 			continue;
// 		// Edge e = make_pair(min(vMin, u), max(vMin, u));
// 		// edge2sup.erase(e);
// 		DEdges.push_back(make_pair(min(u, v), max(u, v)));
// 		set<int> WT;
// 		set_intersection(g->AdjList[v].begin(), g->AdjList[v].end(), g->AdjList[u].begin(), g->AdjList[u].end(), inserter(WT, WT.begin()));
// 		for (int w : WT)
// 		{
// 			if (w < u || removeV.count(w))
// 				continue;

// 			Edge e = make_pair(u, w);
// 			auto it = edge2sup.find(e);
// 			if (it != edge2sup.end())
// 			{

// 				if (it->second < k - 1)
// 				{
// 					return false;
// 				}
// 				else
// 				{
// 					newSupEs.push_back(e);
// 				}
// 			}
// 		}
// 	}
// 	// check whether exists k support edge

// 	for (Edge e : newSupEs)
// 	{
// 		edge2sup[e]--;
// 	}
// 	return true;
// }
bool checkInd(int v, DataGraph *g, set<int> &removeV, int k, unordered_map<int, unordered_map<int, int>> &edge2sup)
{
	vector<Edge> newSupEs;
	sort(g->AdjList[v].begin(), g->AdjList[v].end());
	for (int u : g->AdjList[v])
	{
		if (removeV.count(u))
			continue;
		// Edge e = make_pair(min(vMin, u), max(vMin, u));
		// edge2sup.erase(e);
		sort(g->AdjList[u].begin(), g->AdjList[u].end());
		set<int> WT;
		set_intersection(g->AdjList[v].begin(), g->AdjList[v].end(), g->AdjList[u].begin(), g->AdjList[u].end(), inserter(WT, WT.begin()));
		for (int w : WT)
		{
			if (w < u || removeV.count(w))
				continue;

			Edge e = make_pair(u, w);

			if (edge2sup[u][w] < k - 1)
			{
				return false;
			}
			else
			{
				newSupEs.push_back(e);
			}
		}
	}

	for (Edge e : newSupEs)
	{
		int u = e.first;
		int v = e.second;
		edge2sup[u][v]--;
		edge2sup[v][u]--;
	}
	for (auto e : edge2sup[v])
	{
		int u = e.first;
		edge2sup[u].erase(v);
	}
	edge2sup.erase(v);

	return true;
}
bool checkInd(int v, int k, unordered_map<int, unordered_map<int, int>> &edge2sup)
{
	vector<Edge> newSupEs;
	auto &t = edge2sup[v];
	for (auto &U : t)
	{
		int u = U.first;
		for (auto &t1 : edge2sup[u])
		{
			int w = t1.first;
			if (w < u || !t.count(w))
				continue;

			Edge e = make_pair(u, w);

			if (edge2sup[u][w] < k - 1)
			{
				return false;
			}
			else
			{
				newSupEs.push_back(e);
			}
		}
	}

	for (Edge e : newSupEs)
	{
		int u = e.first;
		int v = e.second;
		edge2sup[u][v]--;
		edge2sup[v][u]--;
	}
	for (auto e : edge2sup[v])
	{
		int u = e.first;
		edge2sup[u].erase(v);
	}
	edge2sup.erase(v);

	return true;
}
// delete edges
// bool checkInd(int v, int k, unordered_map<int, unordered_map<int, int>> &edge2sup)
// {
// 	unordered_map<int, unordered_map<int, int>> edge2supCopy = edge2sup;
// 	unordered_map<int, set<int>> vertex2Nei;
// 	vector<Edge> uwDelCand;
// 	unordered_map<int, int> uwDel;
// 	vector<Edge> mutil1;
// 	unordered_map<int, unordered_map<int, int>> edgeNum;
// 	set<int> vNei;
// 	for (auto e : edge2sup[v])
// 	{
// 		int u = e.first;
// 		vNei.insert(u);
// 	}
// 	for (int u : vNei)
// 	{
// 		set<int> uNei;
// 		for (auto e : edge2sup[u])
// 		{
// 			int w = e.first;
// 			uNei.insert(w);
// 		}
// 		vertex2Nei.emplace(u, uNei);
// 		set<int> comW;
// 		set_intersection(uNei.begin(), uNei.end(), vNei.begin(), vNei.end(), inserter(comW, comW.begin()));
// 		for (int w : comW)
// 		{
// 			if (u > w)
// 				continue;
// 			--edge2supCopy[w][u];
// 			if (--edge2supCopy[u][w] < k - 2)
// 			{
// 				uwDelCand.push_back(make_pair(u, w));
// 			}
// 		}
// 	}
// 	vector<Edge> visit;
// 	while (!uwDelCand.empty())
// 	{
// 		Edge e = uwDelCand.back();
// 		uwDelCand.pop_back();
// 		int u = e.first;
// 		int w = e.second;
// 		if (!edge2supCopy[u].count(w))
// 			continue;
// 		visit.push_back(e);
// 		set<int> wNei;
// 		for (auto e : edge2supCopy[w])
// 		{
// 			int x = e.first;
// 			wNei.insert(x);
// 		}
// 		set<int> uNei;

// 		for (auto e : edge2supCopy[u])
// 		{
// 			int x = e.first;
// 			uNei.insert(x);
// 		}

// 		set<int> comX;
// 		set_intersection(wNei.begin(), wNei.end(), uNei.begin(), uNei.end(), inserter(comX, comX.begin()));
// 		if (comX.size() <= 1) //=1 when just contain v
// 			return false;

// 		int cnt = comX.size(); // 可能v1w1只有一个x且是属于vNei的，但可能在后面viwi删除后，那个点与v1w1构成的两条边里有边<k, 从而v1w1又变得不可删了
// 		for (int x : comX)
// 		{
// 			if (x == v)
// 				continue;

// 			--edge2supCopy[x][u];
// 			if (--edge2supCopy[u][x] < k - 2)
// 			{
// 				uwDelCand.push_back(make_pair(u, w));
// 			}
// 			--edge2supCopy[x][w];
// 			if (--edge2supCopy[w][x] < k - 2)
// 			{
// 				uwDelCand.push_back(make_pair(u, w));
// 			}
// 		}

// 		edge2supCopy[u].erase(w);
// 		edge2supCopy[w].erase(u);
// 	}
// 	for (int u : vNei)
// 	{
// 		edge2supCopy[u].erase(v);
// 	}
// 	edge2supCopy.erase(v);
// 	while (!visit.empty())
// 	{
// 		Edge e = uwDelCand.back();
// 		uwDelCand.pop_back();
// 		int u = e.first;
// 		int w = e.second;
// 		set<int> wNei;
// 		for (auto e : edge2supCopy[w])
// 		{
// 			int x = e.first;
// 			wNei.insert(x);
// 		}
// 		set<int> uNei;

// 		for (auto e : edge2supCopy[u])
// 		{
// 			int x = e.first;
// 			uNei.insert(x);
// 		}

// 		set<int> comX;
// 		set_intersection(wNei.begin(), wNei.end(), uNei.begin(), uNei.end(), inserter(comX, comX.begin()));
// 		if (comX.empty()) //=1 when just contain v
// 			return false;
// 	}

// 	edge2sup = edge2supCopy;

// 	return true;
// }

// bool checkInd(int v, int k, unordered_map<int, unordered_map<int, int>> &edge2sup)
// {
// 	// unordered_map<int, unordered_map<int, int>> edge2supCopy=edge2sup;
// 	vector<Edge> uwDel;
// 	vector<Edge> mutil1;
// 	unordered_map<int, unordered_map<int, int>> edgeNum;
// 	set<int> vNei;
// 	for (auto e : edge2sup[v])
// 	{
// 		int u = e.first;
// 		vNei.insert(u);
// 	}
// 	for (auto e : edge2sup[v])
// 	{

// 		int u = e.first;
// 		set<int> uNei;
// 		for (auto e : edge2sup[u])
// 		{
// 			int w = e.first;
// 			uNei.insert(w);
// 		}
// 		set<int> comW;
// 		set_intersection(uNei.begin(), uNei.end(), vNei.begin(), vNei.end(), inserter(comW, comW.begin()));
// 		for (int w : comW)
// 		{
// 			if (u > w)
// 				continue;
// 			if (edge2sup[u][w] < k - 1)
// 			{
// 				set<int> wNei;
// 				for (auto e : edge2sup[w])
// 				{
// 					int x = e.first;
// 					wNei.insert(x);
// 				}
// 				set<int> comX;
// 				set_intersection(wNei.begin(), wNei.end(), uNei.begin(), uNei.end(), inserter(comX, comX.begin()));
// 				if (comX.size() <= 1)
// 					return false;

// 				for (int x : comX)
// 				{
// 					if (x == v)
// 						continue;
// 					// if (edge2sup[u][x] < k - 1 || edge2sup[w][x] < k - 1)
// 					// {
// 					// 	return false;
// 					// }
// 					int sup = edge2sup[u][x];
// 					int i = min(u, x);
// 					int a = max(u, x);
// 					auto t = edgeNum.find(i);
// 					if (t == edgeNum.end()) // not ux
// 						edgeNum.emplace(i, unordered_map<int, int>{{a, sup}});
// 					else
// 					{
// 						if (!t->second.count(a)) // have u but no x
// 							edgeNum[i].emplace(a, sup);
// 					}
// 					if (--edgeNum[i][a] < k - 2)
// 					{
// 						return false;
// 					}
// 					sup = edge2sup[w][x];
// 					i = min(w, x);
// 					a = max(w, x);
// 					t = edgeNum.find(i);
// 					if (t == edgeNum.end())
// 						edgeNum.emplace(i, unordered_map<int, int>{{a, sup}});
// 					else
// 					{
// 						if (!t->second.count(a))
// 							edgeNum[i].emplace(a, sup);
// 					}

// 					if (--edgeNum[i][a] < k - 2)
// 					{
// 						return false;
// 					}
// 				}
// 				uwDel.push_back(make_pair(u, w));
// 			}
// 			else
// 				mutil1.push_back(make_pair(u, w));
// 		}
// 	}
// 	for (auto e : edge2sup[v])
// 	{
// 		int u = e.first;
// 		edge2sup[u].erase(v);
// 	}
// 	edge2sup.erase(v);

// 	for (Edge e : uwDel)
// 	{
// 		int u = e.first;
// 		int v = e.second;
// 		if(!edge2sup.count(u)||!edge2sup[u].count(v)) cout << "edge not exist" << endl;
// 		edge2sup[u].erase(v);
// 		edge2sup[v].erase(u);
// 	}
// 	for (Edge e : mutil1)
// 	{
// 		int u = e.first;
// 		int v = e.second;
// 		if(!edge2sup.count(u)||!edge2sup[u].count(v)) cout << "edge not exist" << endl;
// 		--edge2sup[u][v];
// 		edge2sup[v][u]--;
// 	}
// 	for (auto t : edgeNum)
// 	{
// 		int u = t.first;
// 		for (auto e : t.second)
// 		{
// 			if(!edge2sup.count(u)||!edge2sup[u].count(e.first)) cout << "edge not exist" << endl;
// 			edge2sup[u][e.first] = e.second;
// 			edge2sup[e.first][u] = e.second;
// 		}
// 	}

// 	return true;
// }
bool checkInd(int v, DataGraph *g, set<int> &removeV, int k, map<Edge, int> &edge2sup, vector<Edge> &kSupEs)
{
	vector<Edge> DEdges;
	vector<Edge> newKSupEs;
	vector<Edge> newSupEs;
	sort(g->AdjList[v].begin(), g->AdjList[v].end());
	for (int u : g->AdjList[v])
	{
		if (removeV.count(u))
			continue;
		// Edge e = make_pair(min(vMin, u), max(vMin, u));
		// edge2sup.erase(e);
		DEdges.push_back(make_pair(min(u, v), max(u, v)));
		sort(g->AdjList[u].begin(), g->AdjList[u].end());
		set<int> WT;
		set_intersection(g->AdjList[v].begin(), g->AdjList[v].end(), g->AdjList[u].begin(), g->AdjList[u].end(), inserter(WT, WT.begin()));
		for (int w : WT)
		{
			if (w < u || removeV.count(w))
				continue;

			Edge e = make_pair(u, w);
			auto it = edge2sup.find(e);
			if (it != edge2sup.end())
			{

				if (it->second < k - 1)
				{
					return false;
				}
				else
				{
					if (it->second == k - 1)
						newKSupEs.push_back(e);
					else
						newSupEs.push_back(e);
				}
			}
		}
	}
	// check whether exists k support edge
	int cnt = 0;
	if (newKSupEs.empty() && DEdges.size() >= kSupEs.size())
	{
		for (auto &e : kSupEs)
		{
			if (!count(DEdges.begin(), DEdges.end(), e))
			{
				cnt++;
				break;
			}
		}
		if (cnt == 0)
			return false;
	}

	for (auto &e : DEdges)
	{
		auto t = find(kSupEs.begin(), kSupEs.end(), e);
		if (t != kSupEs.end())
		{
			kSupEs.erase(t);
		}
	}
	for (Edge e : newKSupEs)
	{
		kSupEs.push_back(e);
		edge2sup[e]--;
	}
	for (Edge e : newSupEs)
	{
		edge2sup[e]--;
	}
	return true;
}
bool checkInd(int v, int k, unordered_map<int, unordered_map<int, int>> &edge2sup, vector<Edge> &kSupEs)
{
	vector<Edge> DEdges;
	vector<Edge> newKSupEs;
	vector<Edge> newSupEs;
	// sort(g->AdjList[v].begin(), g->AdjList[v].end());
	set<int> vNei;
	for (auto &u : edge2sup[v])
	{
		vNei.insert(u.first);
	}

	for (int u : vNei)
	{
		// Edge e = make_pair(min(vMin, u), max(vMin, u));
		// edge2sup.erase(e);
		DEdges.push_back(make_pair(min(u, v), max(u, v)));
		set<int> uNei;
		for (auto &w : edge2sup[u])
		{
			uNei.insert(w.first);
		}
		// sort(g->AdjList[u].begin(), g->AdjList[u].end());
		set<int> WT;
		set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), inserter(WT, WT.begin()));
		for (int w : WT)
		{
			if (w < u)
				continue;

			Edge e = make_pair(u, w);

			if (edge2sup[u][w] < k - 1)
			{
				return false;
			}
			else
			{
				if (edge2sup[u][w] == k - 1)
					newKSupEs.push_back(e);
				else
					newSupEs.push_back(e);
			}
		}
	}
	// check whether exists k support edge
	int cnt = 0;
	if (newKSupEs.empty() && DEdges.size() >= kSupEs.size())
	{
		for (auto &e : kSupEs)
		{
			if (!count(DEdges.begin(), DEdges.end(), e))
			{
				cnt++;
				break;
			}
		}
		if (cnt == 0)
			return false;
	}

	for (auto &e : DEdges)
	{
		auto t = find(kSupEs.begin(), kSupEs.end(), e);
		if (t != kSupEs.end())
		{
			kSupEs.erase(t);
		}
	}
	for (Edge e : newKSupEs)
	{
		kSupEs.push_back(e);
		int u = e.first, w = e.second;
		edge2sup[u][w]--;
		edge2sup[w][u]--;
	}
	for (Edge e : newSupEs)
	{
		int u = e.first, w = e.second;
		edge2sup[u][w]--;
		edge2sup[w][u]--;
	}
	for (auto &e : edge2sup[v])
	{
		edge2sup[e.first].erase(v);
	}
	edge2sup.erase(v);
	return true;
}
// 定义支配关系
bool dominates(const SkyGroupCand &a, const SkyGroupCand &b)
{
	// a 支配 b 的条件
	return (a.k >= b.k && a.rel >= b.rel && a.rho >= b.rho) &&
		   (a.k > b.k || a.rel > b.rel || a.rho > b.rho);
}

// 计算 Skyline
void skyline(vector<SkyGroupCand> &result, const vector<SkyGroupCand> &points)
{
	int n = result.size();

	for (const auto &point : points)
	{
		bool notdom = true;
		for (int i = 0; i < n; i++)
		{
			if (dominates(result[i], point))
			{
				notdom = false;
				break;
			}
		}
		if (notdom == true)
		{
			result.push_back(point);
			// cout << "rel: " << point.rel << " k: " << point.k << " rho: " << point.rho << " maxRV: ";
			// for (auto v : point.vertices)
			// 	cout << v << ", ";
			// cout << endl;
		}
	}
}
bool firstMaxRhoLesser(unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho)
{

	queue<int> lessMRhoV;
	unordered_set<int> removed; // 记录已删除的点
	for (auto &p : vert2weightCopy)
	{
		if (p.second <= maxRho)
		{
			lessMRhoV.push(p.first);
			removed.insert(p.first);
		}
	}

	while (!lessMRhoV.empty())
	{
		int v = lessMRhoV.front();
		lessMRhoV.pop();
		for (auto &e : simG[v])
		{
			if (vert2weightCopy.count(e.first))
			{
				vert2weightCopy[e.first] -= (e.second + denStyle);
				if (vert2weightCopy[e.first] <= maxRho && !removed.count(e.first))
				{
					lessMRhoV.push(e.first);
					removed.insert(e.first);
				}
			}
		}
		vert2weightCopy.erase(v);
	}
	if (vert2weightCopy.empty())
		return true;
	else
		return false;
}
void dfs(int node, unordered_set<int> &visited, unordered_map<int, unordered_map<int, int>> &G)
{
	visited.insert(node);

	for (auto &e : G[node])
	{
		int neighbor = e.first;
		if (visited.find(neighbor) == visited.end())
		{
			dfs(neighbor, visited, G);
		}
	}
}
void firstMaxRhoLesser(unordered_map<int, unordered_map<int, int>> &shrinkG, int k, int ccId, DataGraph *G, unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho)
{
	queue<int> lessMRhoV;
	unordered_set<int> removed; // 记录已删除的点
	for (auto &p : vert2weightCopy)
	{
		if (p.second <= maxRho)
		{
			lessMRhoV.push(p.first);
			removed.insert(p.first);
		}
	}
	if(removed.empty()) return;
	while (!lessMRhoV.empty())
	{
		int v = lessMRhoV.front();
		lessMRhoV.pop();
		for (auto &e : simG[v])
		{
			if (vert2weightCopy.count(e.first))
			{
				vert2weightCopy[e.first] -= (e.second + denStyle);
				if (vert2weightCopy[e.first] <= maxRho && !removed.count(e.first))
				{
					lessMRhoV.push(e.first);
					removed.insert(e.first);
				}
			}
		}
		vert2weightCopy.erase(v);
	}
	if(vert2weightCopy.empty()) return;
	// is k-truss?
	unordered_map<int, unordered_map<int, int>> e2sup;
	queue<Edge> supLesserKE;
	G->support(k, ccId, e2sup, supLesserKE, removed);
	while (!supLesserKE.empty())
	{
		Edge e = supLesserKE.front();
		supLesserKE.pop();
		int u = e.first, v = e.second;

		e2sup[u].erase(v); // remove edge
		e2sup[v].erase(u);
		for (auto w2sup : e2sup[u])
		{
			int w = w2sup.first;
			if (e2sup[v].count(w))
			{
				if (e2sup[u][w] >= k - 2) // if < k-2, which means it already in supLesserKE
				{
					Edge e1 = make_pair(min(w, u), max(w, u));
					--e2sup[u][w];
					if (--e2sup[w][u] < k - 2)
					{
						supLesserKE.push(e1);
					}
				}
				if (e2sup[v][w] >= k - 2)
				{
					Edge e2 = make_pair(min(w, v), max(w, v));
					--e2sup[v][w];
					if (--e2sup[w][v] < k - 2)
					{
						supLesserKE.push(e2);
					}
				}
			}
		}
	}
	if(e2sup.empty()) return;// disconnected
	// connected?
	unordered_set<int> visited;
	int startV = e2sup.begin()->first;
	dfs(startV, visited, e2sup);
	int minSup = INT32_MAX;
	if (visited.size() == e2sup.size())
	{
		// minSup = k-2?
		for (auto &U : e2sup)
		{
			int u = U.first;
			for (auto &V : U.second)
			{
				int v = V.first;
				if (u < v && V.second < minSup)
				{
					minSup = V.second;
				}
			}
		}
		// cout << "kshrink" << endl;
	}
	if (minSup == k - 2)
	{
		shrinkG = e2sup;
		cout << "shrink" << endl;
	}
}
// void firstMaxRhoLesser(unordered_map<int, unordered_map<int, int>> &shrinkG, int k, int ccId, DataGraph *G, unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho)
// {

// 	unordered_set<int> removed; // 记录已删除的点
// 	for (auto &p : vert2weightCopy)
// 	{
// 		if (p.second < maxRho)
// 		{
// 			removed.insert(p.first);
// 		}
// 	}
// 	// is k-truss?
// 	unordered_map<int, unordered_map<int, int>> e2sup;
// 	queue<Edge> supLesserKE;
// 	G->support(k, ccId, e2sup, supLesserKE, removed);
// 	while (!supLesserKE.empty())
// 	{
// 		Edge e = supLesserKE.front();
// 		supLesserKE.pop();
// 		int u = e.first, v = e.second;

// 		e2sup[u].erase(v); // remove edge
// 		e2sup[v].erase(u);
// 		for (auto w2sup : e2sup[u])
// 		{
// 			int w = w2sup.first;
// 			if (e2sup[v].count(w))
// 			{
// 				if (e2sup[u][w] >= k - 2) // if < k-2, which means it already in supLesserKE
// 				{
// 					Edge e1 = make_pair(min(w, u), max(w, u));
// 					--e2sup[u][w];
// 					if (--e2sup[w][u] < k - 2)
// 					{
// 						supLesserKE.push(e1);
// 					}
// 				}
// 				if (e2sup[v][w] >= k - 2)
// 				{
// 					Edge e2 = make_pair(min(w, v), max(w, v));
// 					--e2sup[v][w];
// 					if (--e2sup[w][v] < k - 2)
// 					{
// 						supLesserKE.push(e2);
// 					}
// 				}
// 			}
// 		}
// 	}
// 	// connected?
// 	unordered_set<int> visited;
// 	int startV = e2sup.begin()->first;
// 	dfs(startV, visited, e2sup);
// 	int minSup = INT32_MAX;
// 	if (visited.size() == e2sup.size())
// 	{
// 		// minSup = k-2?
// 		for (auto &U : e2sup)
// 		{
// 			int u = U.first;
// 			for (auto &V : U.second)
// 			{
// 				int v = V.first;
// 				if (u < v && V.second < minSup)
// 				{
// 					minSup = V.second;
// 				}
// 			}
// 		}
// 		cout << "kshrink" << endl;
// 	}
// 	if (minSup == k - 2)
// 	{
// 		shrinkG = e2sup;
// 		cout << "shrink" << endl;
// 	}
// 	while (!removed.empty())
// 	{
// 		visited = removed;
// 		removed.clear();
// 		for (int v : visited)
// 		{
// 			for (auto &e : simG[v])
// 			{
// 				if (vert2weightCopy.count(e.first))
// 				{
// 					vert2weightCopy[e.first] -= (e.second + denStyle);
// 					if (vert2weightCopy[e.first] < maxRho && !visited.count(e.first))
// 					{
// 						removed.insert(e.first);
// 					}
// 				}
// 			}
// 			vert2weightCopy.erase(v);
// 		}
// 		// update support
// 		for (int v : removed)
// 		{
// 			for (auto U : e2sup[v])
// 			{
// 				int u = U.first;
// 				for (auto &W : e2sup[u])
// 				{
// 					int w = W.first;
// 					if (w < u)
// 						continue;
// 					if (e2sup[v].count(w))
// 					{
// 						--e2sup[u][w];
// 						if (--e2sup[w][u] < k - 2)
// 						{
// 							Edge e = make_pair(u, w);
// 							supLesserKE.push(e);
// 						}
// 					}
// 				}
// 				e2sup[u].erase(v);
// 			}
// 			e2sup.erase(v);
// 		}
// 		while (!supLesserKE.empty())
// 		{
// 			Edge e = supLesserKE.front();
// 			supLesserKE.pop();
// 			int u = e.first, v = e.second;

// 			e2sup[u].erase(v); // remove edge
// 			e2sup[v].erase(u);
// 			for (auto w2sup : e2sup[u])
// 			{
// 				int w = w2sup.first;
// 				if (e2sup[v].count(w))
// 				{
// 					if (e2sup[u][w] >= k - 2) // if < k-2, which means it already in supLesserKE
// 					{
// 						Edge e1 = make_pair(min(w, u), max(w, u));
// 						--e2sup[u][w];
// 						if (--e2sup[w][u] < k - 2)
// 						{
// 							supLesserKE.push(e1);
// 						}
// 					}
// 					if (e2sup[v][w] >= k - 2)
// 					{
// 						Edge e2 = make_pair(min(w, v), max(w, v));
// 						--e2sup[v][w];
// 						if (--e2sup[w][v] < k - 2)
// 						{
// 							supLesserKE.push(e2);
// 						}
// 					}
// 				}
// 			}
// 		}
// 		// connected?
// 		unordered_set<int> visited;
// 		int startV = e2sup.begin()->first;
// 		dfs(startV, visited, e2sup);
// 		int minSup = INT32_MAX;
// 		if (visited.size() == e2sup.size())
// 		{
// 			// minSup = k-2?
// 			for (auto &U : e2sup)
// 			{
// 				int u = U.first;
// 				for (auto &V : U.second)
// 				{
// 					int v = V.first;
// 					if (u < v && V.second < minSup)
// 					{
// 						minSup = V.second;
// 					}
// 				}
// 			}
// 			cout << "kshrink" << endl;
// 		}
// 		if (minSup == k - 2)
// 		{
// 			shrinkG = e2sup;
// 			cout << "shrink" << endl;
// 		}
// 	}
// }
bool secondMaxRhoLesser(unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho, int delV)
{

	queue<int> lessMRhoV;
	unordered_set<int> removed;
	if (vert2weightCopy.count(delV))
	{
		for (auto &e : simG[delV])
		{
			if (vert2weightCopy.count(e.first))
			{
				vert2weightCopy[e.first] -= (e.second + denStyle);
				if (vert2weightCopy[e.first] <= maxRho)
				{
					lessMRhoV.push(e.first);
					removed.insert(e.first);
				}
			}
		}
		vert2weightCopy.erase(delV);
	}
	while (!lessMRhoV.empty())
	{
		int v = lessMRhoV.front();
		lessMRhoV.pop();
		for (auto &e : simG[v])
		{
			if (vert2weightCopy.count(e.first))
			{
				vert2weightCopy[e.first] -= (e.second + denStyle);
				if (vert2weightCopy[e.first] <= maxRho && !removed.count(e.first))
				{
					lessMRhoV.push(e.first);
					removed.insert(e.first);
				}
			}
		}
		vert2weightCopy.erase(v);
	}
	if (vert2weightCopy.empty())
		return true;
	else
		return false;
}
void updateBcc(set<int> &cv, unordered_map<int, vector<int>> &bccId2verts, set<int> &cv1, unordered_map<int, vector<int>> &bccId2verts1, unordered_map<int, vector<int>> &vert2bccId, int bccId, int &bccIndex)
{
	if (bccId2verts1.size() == 0) // if the bcc is empty after deleting a non-cut which means that the whole graph is deleted be empty
	{
		// if (bccId2verts.size() != 1)
		// bccId2verts[bccId].clear();
		bccId2verts.erase(bccId);
	}
	else if (bccId2verts1.size() == 1) // means that this bcc does not develop any cut vertices after deleting a non-cut
	{
		auto t = bccId2verts1.begin()->second;
		if (t.size() == 1) // the only vertex must is cut vertex in the whole graph
		{
			int v = t[0];
			auto tp = find(vert2bccId[v].begin(), vert2bccId[v].end(), bccId);
			vert2bccId[v].erase(tp);
			if (vert2bccId[v].size() < 2)
				cv.erase(v);
			if (bccId2verts.size() > 1)
				// bccId2verts[bccId].clear();
				bccId2verts.erase(bccId);
		}
	}
	else
	{
		int maxSizeId = 0;
		for (auto t : bccId2verts1)
		{
			if (t.second.size() > bccId2verts[maxSizeId].size())
				maxSizeId = t.first;
		}

		auto &remainBcc = bccId2verts1[maxSizeId];
		unordered_set<int> cvInRBcc;
		for (int v : cv1)
		{
			if (count(remainBcc.begin(), remainBcc.end(), v))
				cvInRBcc.insert(v);
		}
		bccId2verts1.erase(maxSizeId);
		unordered_map<int, vector<int>> vert2bccId1;
		auto &vs = bccId2verts[bccId];
		for (auto t : bccId2verts1)
		{

			// 0-1,0-5,1-2,1-5,2-3,3-4,3-5,4-5,5-6,5-7,6-7;cut 5,bcc1:0,1,2,3,4,5;bcc2:5,6,7;delete 2;bcc1:0,1,5;bcc2:3,4,5;bcc3:5,6,7
			bccId2verts.emplace(bccIndex, t.second);
			for (auto v : t.second)
			{

				vector<int> &oldBccId = vert2bccId[v];
				if (cv1.count(v))
				{
					if (!cvInRBcc.count(v))
					{
						auto tmp = find(vs.begin(), vs.end(), v);
						if (tmp != vs.end()) // if ==end, which means the cut is a new
						{
							vs.erase(tmp);
						}
					}
					auto tmp = vert2bccId1.find(v);
					if (tmp == vert2bccId1.end())
					{
						vert2bccId1.emplace(v, vector<int>{bccIndex});
					}
					else
						vert2bccId1[v].push_back(bccIndex);
				}
				else // the new non-cut only belongs to a new bcc, so can immediate update
				{
					// the common vertex with the maxmimum new bcc must be the new cut vertex, so the new non-cut vertices can immediate delete from the maxmimum new bcc

					auto tmp = find(vs.begin(), vs.end(), v);
					vs.erase(tmp);

					auto tp = find(oldBccId.begin(), oldBccId.end(), bccId);
					if (tp != oldBccId.end())
					{
						*tp = bccIndex;
					}
				}
			}
			bccIndex++;
		}
		cv.insert(cv1.begin(), cv1.end());
		for (int v : cv1) // the new cut vertices may from cut and non-cut in the whole graph
		{
			vector<int> &oldBccId = vert2bccId[v];
			if (cvInRBcc.count(v)) // save the origin bccId and add new id
			{
				for (int u : vert2bccId1[v])
					oldBccId.push_back(u);
			}
			else
			{
				if (oldBccId.size() == 1)
				{
					vert2bccId[v] = vert2bccId1[v];
				}
				else
				{
					auto tp = find(oldBccId.begin(), oldBccId.end(), bccId);
					if (tp != oldBccId.end())
					{
						oldBccId.erase(tp);
					}
					for (int v : vert2bccId1[v])
						oldBccId.push_back(v);
				}
			}
		}
	}
}
void BCCUtil(set<int> &cv, int u, unordered_map<int, int> &disc, unordered_map<int, int> &low, vector<int> &st,
			 unordered_map<int, int> &parent, unordered_map<int, vector<int>> &bccs, int &time, int &index, unordered_map<int, unordered_map<int, int>> &g)
{
	disc[u] = low[u] = ++time;
	int children = 0;
	st.push_back(u);
	for (auto e : g[u])
	{
		int v = e.first;
		if (!disc.count(v))
			continue;
		if (disc[v] == -1)
		{
			children++;
			parent[v] = u;

			BCCUtil(cv, v, disc, low, st, parent, bccs, time, index, g);

			low[u] = min(low[u], low[v]);

			if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u]))
			{
				cv.insert(u);
				vector<int> bcc;
				int x;
				do
				{
					x = st.back();
					bcc.push_back(x);
					st.pop_back();
				} while (x != v);
				// while (st.back() != u)
				// {
				//     bcc.push_back(st.back());
				//     st.pop_back();
				// }
				bcc.push_back(u); // Add u (articulation point)

				bccs[index++] = bcc;
			}
		}

		else if (v != parent[u])
		{
			low[u] = min(low[u], disc[v]);
		}
	}
}
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, set<int> &subG, unordered_map<int, unordered_map<int, int>> &g)
{
	unordered_map<int, int> disc;
	unordered_map<int, int> low;
	unordered_map<int, int> parent;
	int time = 0;
	int index = 0;
	// vector<vector<int>> bccs;
	vector<int> st;

	// get neighbors of delete vertex
	//  int dv;
	// vector<int> trueSub;
	// Initialize disc and low, and parent arrays
	for (int i : subG)
	{
		disc[i] = NIL;
		low[i] = NIL;
		parent[i] = NIL;
		// trueSub.push_back(i);
	}
	for (int i : subG)
	{
		auto t = disc.find(i);
		if (t != disc.end() && t->second == NIL)
			BCCUtil(cv, i, disc, low, st, parent, bccs, time, index, g);

		// If stack is not empty, pop all edges from stack
		vector<int> bcc;
		while (st.size() > 0)
		{

			bcc.push_back(st.back());
			st.pop_back();
		}
		if (!bcc.empty())
		{
			bccs[index++] = bcc;
		}
	}
}
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, unordered_map<int, unordered_map<int, int>> &g)
{
	unordered_map<int, int> disc;
	unordered_map<int, int> low;
	unordered_map<int, int> parent;
	int time = 0;
	int index = 0;
	// vector<vector<int>> bccs;
	vector<int> st;

	// get neighbors of delete vertex
	//  int dv;
	// vector<int> trueSub;
	// Initialize disc and low, and parent arrays
	for (int i : subG)
	{
		disc[i] = NIL;
		low[i] = NIL;
		parent[i] = NIL;
		// trueSub.push_back(i);
	}
	for (int i : subG)
	{
		auto t = disc.find(i);
		if (t != disc.end() && t->second == NIL)
			BCCUtil(cv, i, disc, low, st, parent, bccs, time, index, g);

		// If stack is not empty, pop all edges from stack
		vector<int> bcc;
		while (st.size() > 0)
		{

			bcc.push_back(st.back());
			st.pop_back();
		}
		if (!bcc.empty())
		{
			bccs[index++] = bcc;
		}
	}
}
void BCCUtil(set<int> &cv, int u, unordered_map<int, int> &disc, unordered_map<int, int> &low, vector<int> &st,
			 unordered_map<int, int> &parent, unordered_map<int, vector<int>> &bccs, set<int> &delV, int &time, int &index, DataGraph *relSubG)
{
	disc[u] = low[u] = ++time;
	int children = 0;
	st.push_back(u);
	for (int v : relSubG->AdjList[u])
	{
		if (!disc.count(v) || delV.count(v))
			continue;
		if (disc[v] == -1)
		{
			children++;
			parent[v] = u;

			BCCUtil(cv, v, disc, low, st, parent, bccs, delV, time, index, relSubG);

			low[u] = min(low[u], low[v]);

			if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u]))
			{
				cv.insert(u);
				vector<int> bcc;
				int x;
				do
				{
					x = st.back();
					bcc.push_back(x);
					st.pop_back();
				} while (x != v);
				// while (st.back() != u)
				// {
				//     bcc.push_back(st.back());
				//     st.pop_back();
				// }
				bcc.push_back(u); // Add u (articulation point)

				bccs[index++] = bcc;
			}
		}

		else if (v != parent[u])
		{
			low[u] = min(low[u], disc[v]);
		}
	}
}
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, set<int> &delV, int dv, DataGraph *relSubG)
{
	unordered_map<int, int> disc;
	unordered_map<int, int> low;
	unordered_map<int, int> parent;
	int time = 0;
	int index = 0;
	// vector<vector<int>> bccs;
	vector<int> st;

	for (int i : subG)
	{
		if (delV.find(i) != delV.end())
			continue;
		disc[i] = NIL;
		low[i] = NIL;
		parent[i] = NIL;
	}

	for (int i : subG)
	{
		auto t = disc.find(i);
		if (t != disc.end() && t->second == NIL)
			BCCUtil(cv, i, disc, low, st, parent, bccs, delV, time, index, relSubG);

		// If stack is not empty, pop all edges from stack
		vector<int> bcc;
		while (st.size() > 0)
		{

			bcc.push_back(st.back());
			st.pop_back();
		}
		if (!bcc.empty())
		{
			bccs[index++] = bcc;
		}
	}
}
void BCCUtil(set<int> &cv, int u, unordered_map<int, int> &disc, unordered_map<int, int> &low, vector<int> &st,
			 unordered_map<int, int> &parent, unordered_map<int, vector<int>> &bccs, set<int> &delV, int &time, int &index, unordered_map<int, unordered_map<int, int>> &edge2sup)
{
	disc[u] = low[u] = ++time;
	int children = 0;
	st.push_back(u);
	for (auto V : edge2sup[u])
	{
		int v = V.first;
		// if (!disc.count(v) || delV.count(v))
		// 	continue;
		if (!disc.count(v))
			continue;
		// if (delV.count(v)){
		// 	cout<<"delV: "<<v<<endl;
		// 	continue;
		// }

		if (disc[v] == -1)
		{
			children++;
			parent[v] = u;

			BCCUtil(cv, v, disc, low, st, parent, bccs, delV, time, index, edge2sup);

			low[u] = min(low[u], low[v]);

			if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u]))
			{
				cv.insert(u);
				vector<int> bcc;
				int x;
				do
				{
					x = st.back();
					bcc.push_back(x);
					st.pop_back();
				} while (x != v);
				// while (st.back() != u)
				// {
				//     bcc.push_back(st.back());
				//     st.pop_back();
				// }
				bcc.push_back(u); // Add u (articulation point)

				bccs[index++] = bcc;
			}
		}

		else if (v != parent[u])
		{
			low[u] = min(low[u], disc[v]);
		}
	}
}
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, set<int> &delV, unordered_map<int, unordered_map<int, int>> &edge2sup)
{
	unordered_map<int, int> disc;
	unordered_map<int, int> low;
	unordered_map<int, int> parent;
	int time = 0;
	int index = 0;
	// vector<vector<int>> bccs;
	vector<int> st;

	for (int i : subG)
	{
		if (delV.find(i) != delV.end())
		{
			cout << "delV: " << i << endl;
			continue;
		}

		disc[i] = NIL;
		low[i] = NIL;
		parent[i] = NIL;
	}

	for (int i : subG)
	{
		auto t = disc.find(i);
		if (t != disc.end() && t->second == NIL)
			BCCUtil(cv, i, disc, low, st, parent, bccs, delV, time, index, edge2sup);

		// If stack is not empty, pop all edges from stack
		vector<int> bcc;
		while (st.size() > 0)
		{

			bcc.push_back(st.back());
			st.pop_back();
		}
		if (!bcc.empty())
		{
			bccs[index++] = bcc;
		}
	}
}