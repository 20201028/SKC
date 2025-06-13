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
#include <cmath>
// #include "DataGraph.h"
// #include "KTruss.h"
#include "Basic.h"
#include "Improve.h"
#include "Best.h"
// #include "Define.h"
// #include "EquiTree.h"

using namespace std;

unordered_map<int, set<int>> id2att;
unordered_map<int, pair<int, int>> id2loc;
unordered_map<int, set<int>> att2vid;
unordered_map<int, unordered_set<int>> oriG;

// map<int, unordered_map<int, double>> simE;
// unordered_set<int> existV;

// int denStyle = 0; // 1:edge add weight; 0: only weight
int relThd = 0;

int ce = 0, cn = 0;

// 计算连通分量的函数
void findConnectedComponents(unordered_set<int> &subG,
							 unordered_map<int, int> &componentMap)
{
	unordered_set<int> visited;
	visited.reserve(subG.size());
	int componentId = 0;
	unordered_map<int, shared_ptr<UNode>> idUFMap;

	// 遍历图中的每个节点
	for (const auto &node : subG)
	{
		shared_ptr<UNode> unode;
		if (idUFMap.count(node))
			unode = idUFMap[node];
		else
		{
			unode = make_shared<UNode>(node);
			makeSet(unode);
			idUFMap.emplace(node, unode);
		}

		for (int neighbor : oriG[node])
		{
			if (visited.count(neighbor) || !subG.count(neighbor))
				continue;
			shared_ptr<UNode> unode1;
			if (idUFMap.count(neighbor))
				unode1 = idUFMap[neighbor];
			else
			{
				unode1 = make_shared<UNode>(neighbor);
				makeSet(unode1);
				idUFMap.emplace(neighbor, unode1);
			}
			unite1(unode, unode1);
		}
		visited.insert(node);
	}
	for (auto &v : idUFMap)
	{
		int vid = v.first;
		int root = find(v.second)->value;
		componentMap.emplace(vid, root);
	}
}
// int getSimEs(double simThreshold, string &dataset, unordered_map<int, set<int>> &rel2verts, map<int, unordered_map<int, double>> &simE, int iq)
// {
// 	vector<int> allV;
// 	unordered_set<int> subG;
// 	int relmax = 0;
// 	// if (iq == 39)
// 	// {
// 	// 	for (auto iter : rel2verts)
// 	// 	{
// 	// 		if (iter.first > 1)
// 	// 			for (auto v : rel2verts[iter.first])
// 	// 			{
// 	// 				allV.push_back(v);
// 	// 				subG.insert(v);
// 	// 			}
// 	// 	}
// 	// 	relmax = 2;
// 	// }
// 	// for (auto iter : rel2verts)
// 	// {
// 	// 	if (iter.first > relmax)
// 	// 		relmax = iter.first;
// 	// }
// 	// for (auto v : rel2verts[relmax])
// 	// {
// 	// 	allV.push_back(v);
// 	// }
// 	// else
// 	// {
// 	// 	for (auto iter : rel2verts)
// 	// 	{
// 	// 		for (auto v : iter.second)
// 	// 		{
// 	// 			allV.push_back(v);
// 	// 			subG.insert(v);
// 	// 		}
// 	// 	}
// 	// }
// 	for (auto iter : rel2verts)
// 	{
// 		for (auto v : iter.second)
// 		{
// 			allV.push_back(v);
// 			subG.insert(v);
// 		}
// 	}
// 	unordered_map<int, int> componentMap;
// 	findConnectedComponents(subG, componentMap);
// 	sort(allV.begin(), allV.end());
// 	// string AttFileName = "DataGraph/" + dataset + "/" + "simG.txt";
// 	// ifstream ain(AttFileName.c_str());

// 	// if (!ain)
// 	// {

// 	// ofstream outFile(AttFileName.c_str());
// 	for (int i = 0; i < allV.size(); i++)
// 	{
// 		int u = allV[i];
// 		// existV.insert(u);

// 		// unordered_set<int> nei;
// 		// for (int kw : id2att[u])
// 		// {
// 		// 	set<int> tempIntersection;
// 		// 	set_intersection(att2vid[kw].begin(), att2vid[kw].end(),
// 		// 					 allV.begin(), allV.end(),
// 		// 					 inserter(tempIntersection, tempIntersection.begin()));
// 		// 	nei.insert(tempIntersection.begin(), tempIntersection.end());
// 		// }
// 		// for (int v : nei)
// 		// {

// 		for (int j = i + 1; j < allV.size(); j++)
// 		{
// 			int v = allV[j];

// 			// if (v<=u||componentMap[v] != componentMap[u])
// 			if (componentMap[v] != componentMap[u])
// 				continue;
// 			set<int> inter;
// 			set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

// 			// Weight w(inter.size(), unin.size());
// 			// Weight thd(0,5);
// 			double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
// 			double w = in / un;
// 			if (inter.size() > 0 && w > simThreshold)
// 			{

// 				// outFile << min(u, v) << '\t' << max(u, v) << '\t' << w << endl;
// 				unordered_map<int, double> ew1;
// 				ew1[v] = w;
// 				auto ret1 = simE.emplace(u, ew1);
// 				if (!ret1.second)
// 				{
// 					simE[u][v] = w;
// 				}
// 				unordered_map<int, double> ew2;
// 				ew2[u] = w;
// 				auto ret2 = simE.emplace(v, ew2);
// 				if (!ret2.second)
// 				{
// 					simE[v][u] = w;
// 				}
// 			}
// 		}
// 	}
// 	// outFile.close();
// 	cout << "Construct Successful." << endl;
// 	return relmax;
// 	// }
// 	// else
// 	// {
// 	// 	int u, v;
// 	// 	double w;
// 	// 	unordered_set<int> existV;
// 	// 	while (ain >> u >> v >> w)
// 	// 	{

// 	// 		if (!count(allV.begin(), allV.end(), u) || !count(allV.begin(), allV.end(), v))
// 	// 			continue;
// 	// 		existV.insert(u);
// 	// 		existV.insert(v);
// 	// 		unordered_map<int, double> ew1;
// 	// 		ew1[v] = w;
// 	// 		auto ret1 = simE.emplace(u, ew1);
// 	// 		if (!ret1.second)
// 	// 		{
// 	// 			simE[u][v] = w;
// 	// 		}
// 	// 		unordered_map<int, double> ew2;
// 	// 		ew2[u] = w;
// 	// 		auto ret2 = simE.emplace(v, ew2);
// 	// 		if (!ret2.second)
// 	// 		{
// 	// 			simE[v][u] = w;
// 	// 		}
// 	// 	}
// 	// 	ain.close();
// 	// 	ofstream outFile(AttFileName.c_str(), std::ios::app);
// 	// 	for (int u : allV)
// 	// 	{

// 	// 		if (!simE.count(u))
// 	// 		{
// 	// 			simE.emplace(u, unordered_map<int, double>());
// 	// 			for (int v : existV)
// 	// 			{
// 	// 				if (componentMap[v] != componentMap[u])
// 	// 					continue;
// 	// 				set<int> inter;
// 	// 				set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

// 	// 				// Weight w(inter.size(), unin.size());
// 	// 				// Weight thd(0,5);
// 	// 				double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
// 	// 				double w = in / un;
// 	// 				if (inter.size() > 0 && w > simThreshold)
// 	// 				{
// 	// 					outFile << min(u, v) << '\t' << max(u, v) << '\t' << w << endl;
// 	// 					unordered_map<int, double> ew1;
// 	// 					ew1[v] = w;
// 	// 					auto ret1 = simE.emplace(u, ew1);
// 	// 					if (!ret1.second)
// 	// 					{
// 	// 						simE[u][v] = w;
// 	// 					}
// 	// 					unordered_map<int, double> ew2;
// 	// 					ew2[u] = w;
// 	// 					auto ret2 = simE.emplace(v, ew2);
// 	// 					if (!ret2.second)
// 	// 					{
// 	// 						simE[v][u] = w;
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 			existV.insert(u);
// 	// 		}
// 	// 		else
// 	// 		{
// 	// 			for (int v : existV)
// 	// 			{
// 	// 				if (u == v || componentMap[v] != componentMap[u] || simE[u].count(v))
// 	// 					continue;
// 	// 				set<int> inter;
// 	// 				set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

// 	// 				// Weight w(inter.size(), unin.size());
// 	// 				// Weight thd(0,5);
// 	// 				double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
// 	// 				double w = in / un;
// 	// 				if (inter.size() > 0 && w > simThreshold)
// 	// 				{
// 	// 					outFile << min(u, v) << '\t' << max(u, v) << '\t' << w << endl;
// 	// 					unordered_map<int, double> ew1;
// 	// 					ew1[v] = w;
// 	// 					auto ret1 = simE.emplace(u, ew1);
// 	// 					if (!ret1.second)
// 	// 					{
// 	// 						simE[u][v] = w;
// 	// 					}
// 	// 					unordered_map<int, double> ew2;
// 	// 					ew2[u] = w;
// 	// 					auto ret2 = simE.emplace(v, ew2);
// 	// 					if (!ret2.second)
// 	// 					{
// 	// 						simE[v][u] = w;
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}
// 	// 	}
// 	// 	outFile.close();
// 	// 	cout << "Load simG successfully!" << endl;
// 	// }
// 	// return relmax;
// }
// void getSimEs(double simThreshold, string &dataset, unordered_map<int, set<int>> &rel2verts, map<int, unordered_map<int, double>> &simE)
// {
// 	vector<int> allV;
// 	unordered_set<int> subG;
// 	for (auto iter : rel2verts)
// 	{
// 		for (auto v : iter.second)
// 		{
// 			allV.push_back(v);
// 			subG.insert(v);
// 		}
// 	}
// 	unordered_map<int, int> componentMap;
// 	findConnectedComponents(subG, componentMap);
// 	sort(allV.begin(), allV.end());
// 	// string AttFileName = "DataGraph/" + dataset + "/" + "simG.txt";
// 	// ifstream ain(AttFileName.c_str());

// 	// if (!ain)
// 	// {

// 	// ofstream outFile(AttFileName.c_str());
// 	for (int i = 0; i < allV.size(); i++)
// 	{
// 		int u = allV[i];

// 		for (int j = i + 1; j < allV.size(); j++)
// 		{
// 			int v = allV[j];

// 			// if (v<=u||componentMap[v] != componentMap[u])
// 			if (componentMap[v] != componentMap[u])
// 				continue;
// 			set<int> inter;
// 			set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

// 			// Weight w(inter.size(), unin.size());
// 			// Weight thd(0,5);
// 			double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
// 			double w = in / un;
// 			if (inter.size() > 0 && w > simThreshold)
// 			{

// 				// outFile << min(u, v) << '\t' << max(u, v) << '\t' << w << endl;
// 				unordered_map<int, double> ew1;
// 				ew1[v] = w;
// 				auto ret1 = simE.emplace(u, ew1);
// 				if (!ret1.second)
// 				{
// 					simE[u][v] = w;
// 				}
// 				unordered_map<int, double> ew2;
// 				ew2[u] = w;
// 				auto ret2 = simE.emplace(v, ew2);
// 				if (!ret2.second)
// 				{
// 					simE[v][u] = w;
// 				}
// 			}
// 		}
// 	}
// 	// outFile.close();
// 	cout << "Construct Successful." << endl;
// }
void getMaxRelSimEs(int relmax, string &dataset, unordered_map<int, set<int>> &rel2verts, map<int, unordered_map<int, double>> &simE, string type)
{
	vector<int> allV;
	unordered_set<int> subG;

	for (auto iter : rel2verts)
	{
		if (iter.first >= relmax)
			for (auto v : iter.second)
			{
				allV.push_back(v);
				subG.insert(v);
			}
	}
	unordered_map<int, int> componentMap;
	findConnectedComponents(subG, componentMap);
	sort(allV.begin(), allV.end());

	for (int i = 0; i < allV.size(); i++)
	{
		int u = allV[i];

		for (int j = i + 1; j < allV.size(); j++)
		{
			int v = allV[j];

			if (componentMap[v] != componentMap[u])
				continue;

			double w = 0;
			if (type == "i")
			{
				if (id2att.count(u) && id2att.count(v))
				{
					set<int> inter;
					set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

					double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
					w = in / un;
				}
			}
			else
			{
				if (id2loc.count(u) && id2loc.count(v))
					// w = sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2));
					// w = exp(-pow(id2loc[u].first - id2loc[v].first, 2) - pow(id2loc[u].second - id2loc[v].second, 2));
					w = exp(-0.01 * (sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2))));
			}

			if (w > 0)
			{

				unordered_map<int, double> ew1;
				ew1[v] = w;
				auto ret1 = simE.emplace(u, ew1);
				if (!ret1.second)
				{
					simE[u][v] = w;
				}
				unordered_map<int, double> ew2;
				ew2[u] = w;
				auto ret2 = simE.emplace(v, ew2);
				if (!ret2.second)
				{
					simE[v][u] = w;
				}
			}
		}
	}
	// outFile.close();
	cout << "Construct Successful." << endl;
}
void getSimEs(string &dataset, unordered_map<int, set<int>> &rel2verts, map<int, unordered_map<int, double>> &simE, string type, double simThreshold)
{
	vector<int> allV;
	unordered_set<int> subG;

	for (auto iter : rel2verts)
	{
		for (auto v : iter.second)
		{
			allV.push_back(v);
			subG.insert(v);
		}
	}
	unordered_map<int, int> componentMap;
	findConnectedComponents(subG, componentMap);
	sort(allV.begin(), allV.end());

	for (int i = 0; i < allV.size(); i++)
	{
		int u = allV[i];

		for (int j = i + 1; j < allV.size(); j++)
		{
			int v = allV[j];

			if (componentMap[v] != componentMap[u])
				continue;

			double w = 0;
			if (type == "i")
			{
				if (id2att.count(u) && id2att.count(v))
				{
					set<int> inter;
					set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

					double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
					w = in / un;
				}
			}
			else
			{
				if (id2loc.count(u) && id2loc.count(v))
					// w = sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2));
					w = exp(-0.01 * (sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2))));
			}

			if (w > simThreshold)
			{

				unordered_map<int, double> ew1;
				ew1[v] = w;
				auto ret1 = simE.emplace(u, ew1);
				if (!ret1.second)
				{
					simE[u][v] = w;
				}
				unordered_map<int, double> ew2;
				ew2[u] = w;
				auto ret2 = simE.emplace(v, ew2);
				if (!ret2.second)
				{
					simE[v][u] = w;
				}
			}
		}
	}
	// outFile.close();
	cout << "Construct Successful." << endl;
}
void getKWSubG(string dataset, set<int> &qw, map<int, vector<Edge>, greater<int>> &rel2edge, unordered_map<int, set<int>> &rel2verts)
{
	unordered_map<int, int> vert2rel;
	unordered_set<int> visited;
	for (int w : qw)
	{
		for (int ver : att2vid[w])
		{
			if (!oriG.count(ver) || visited.count(ver))
				continue;
			visited.insert(ver);
			set<int> &att = id2att[ver];
			set<int> inter;
			set_intersection(qw.begin(), qw.end(), att.begin(), att.end(), inserter(inter, inter.begin()));
			int rel = inter.size();
			// if(rel>1) {
			// 	for(int qv : inter){
			// 		cout << qv << endl;
			// 	}
			// }

			if (rel > relThd)
			{
				auto t = rel2verts.find(rel);
				if (t == rel2verts.end())
				{
					rel2verts.emplace(rel, set<int>{ver});
				}
				else
					t->second.insert(ver);

				vert2rel.emplace(ver, rel);
			}
		}
	}
	for (auto &it : vert2rel)
	{
		int src = it.first;
		for (int dst : oriG[src])
			if (src < dst && vert2rel.count(dst))
			{

				Edge e = make_pair(src, dst); // cout<<src<<" "<<dst<<endl;
				int sr = vert2rel[src];
				int dr = vert2rel[dst];
				if (sr > dr)
					swap(sr, dr);
				auto t = rel2edge.find(sr);
				if (t == rel2edge.end())
				{
					vector<Edge> a = {e};
					rel2edge.emplace(sr, a);
				}
				else
					t->second.push_back(e);
			}
	}
}
// void LoadData(int k, string dataset)
// {

// 	string StructFileName = "DataGraph/" + dataset + "/" + "graph.txt";
// 	ifstream sin(StructFileName.c_str());

// 	if (!sin)
// 	{
// 		cout << "Fail to read " << StructFileName << "." << endl;
// 		return;
// 	}
// 	string sline;
// 	while (getline(sin, sline)) // 默认数据集中边不重复
// 	{
// 		if (sline.find('#') != string::npos)
// 			continue;
// 		string src_s = sline.substr(0, sline.find("\t"));
// 		string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
// 		int src = stof(src_s); // string2float
// 		int dst = stof(dst_s);
// 		if (src != dst)
// 		{
// 			auto it = oriG.find(src);
// 			if (it == oriG.end())
// 			{
// 				oriG.emplace(src, unordered_set<int>{dst});
// 			}
// 			else
// 			{
// 				oriG[src].insert(dst);
// 			}
// 			it = oriG.find(dst);
// 			if (it == oriG.end())
// 			{
// 				oriG.emplace(dst, unordered_set<int>{src});
// 			}
// 			else
// 			{
// 				oriG[dst].insert(src);
// 			}
// 		}
// 	}
// 	sin.close();
// 	int d = 0;
// 	for (auto it : oriG)
// 	{
// 		d += it.second.size();
// 	}
// 	cout << "avg degree: " << d / oriG.size() << endl;

// 	unordered_set<int> visited;
// 	unordered_set<int> remove;
// 	for (auto it = oriG.begin(); it != oriG.end();)
// 	{

// 		if (it->second.size() < k - 1)
// 		{
// 			for (int v : it->second)
// 			{
// 				oriG[v].erase(it->first);
// 				if (visited.count(v) && oriG[v].size() < k - 1)
// 				{
// 					remove.insert(v);
// 				}
// 			}
// 			it = oriG.erase(it);
// 		}
// 		else
// 		{
// 			visited.insert(it->first);
// 			++it;
// 		}
// 	}
// 	queue<int> q;
// 	for (int v : remove)
// 	{
// 		q.push(v);
// 	}
// 	while (!q.empty())
// 	{
// 		int v = q.front();
// 		q.pop();

// 		for (int u : oriG[v])
// 		{
// 			oriG[u].erase(v);
// 			if (!remove.count(u) && oriG[u].size() < k - 1)
// 			{
// 				remove.insert(u);
// 				q.push(u);
// 			}
// 		}
// 		oriG.erase(v);
// 	}
// 	string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
// 	ifstream ain(AttFileName.c_str());

// 	if (!ain)
// 	{
// 		cout << "Fail to read " << AttFileName << "." << endl;
// 		return;
// 	}

// 	string aline;
// 	// unordered_map<int, set<int>> nodeAtt;

// 	while (getline(ain, aline))
// 	{
// 		string ver_s = aline.substr(0, aline.find("\t"));
// 		string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
// 		int ver = stof(ver_s); // string2float
// 		if (!oriG.count(ver))
// 			continue;
// 		int start = 0, end = 0;
// 		set<int> att;
// 		int w = 0;
// 		while ((end = att_s.find(",", start)) != string::npos)
// 		{
// 			w = stof(att_s.substr(start, end - start));
// 			att.insert(w);
// 			auto it = att2vid.find(w);
// 			if (it == att2vid.end())
// 			{
// 				att2vid.emplace(w, set<int>{ver});
// 			}
// 			else
// 			{
// 				att2vid[w].insert(ver);
// 			}
// 			start = end + 1;
// 		}
// 		w = stof(att_s.substr(start));
// 		att.insert(w);
// 		auto it = att2vid.find(w);
// 		if (it == att2vid.end())
// 		{
// 			att2vid.emplace(w, set<int>{ver});
// 		}
// 		else
// 		{
// 			att2vid[w].insert(ver);
// 		}

// 		id2att.emplace(ver, att);
// 	}
// 	ain.close();
// 	cout << "Loaded attribute successfully!" << endl;
// 	// findConnectedComponents(oriG, componentMap);
// 	string LocationFileName = "DataGraph/" + dataset + "/" + "user_location.txt";
// 	ifstream lin(LocationFileName.c_str());

// 	if (!lin)
// 	{
// 		cout << "No location information." << endl;
// 		return;
// 	}
// 	string lline;
// 	while (getline(lin, lline))
// 	{
// 		if (lline.find('#') != string::npos)
// 			continue;
// 		size_t first_tab = lline.find("\t");
// 		string user_id = lline.substr(0, first_tab);
// 		size_t second_tab = lline.find("\t", first_tab + 1);
// 		string locs1 = lline.substr(first_tab + 1, second_tab - first_tab - 1);
// 		string locs2 = lline.substr(second_tab + 1);
// 		int uid = stof(user_id); // string2float
// 		if (!oriG.count(uid))
// 			continue;
// 		int loc1 = stof(locs1);
// 		int loc2 = stof(locs2);
// 		id2loc.emplace(uid, make_pair(loc1, loc2));
// 	}
// 	cout << "Loaded dataset successfully!" << endl;
// }
void LoadData(string dataset)
{

	string StructFileName = "DataGraph/" + dataset + "/" + "graph.txt";
	ifstream sin(StructFileName.c_str());

	if (!sin)
	{
		cout << "Fail to read " << StructFileName << "." << endl;
		return;
	}
	string sline;
	while (getline(sin, sline)) // 默认数据集中边不重复
	{
		if (sline.find('#') != string::npos)
			continue;
		string src_s = sline.substr(0, sline.find("\t"));
		string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
		int src = stof(src_s); // string2float
		int dst = stof(dst_s);
		if (src != dst)
		{
			auto it = oriG.find(src);
			if (it == oriG.end())
			{
				oriG.emplace(src, unordered_set<int>{dst});
			}
			else
			{
				oriG[src].insert(dst);
			}
			it = oriG.find(dst);
			if (it == oriG.end())
			{
				oriG.emplace(dst, unordered_set<int>{src});
			}
			else
			{
				oriG[dst].insert(src);
			}
		}
	}
	sin.close();
	int d = 0;
	for (auto it : oriG)
	{
		d += it.second.size();
	}
	cout << "avg degree: " << d / oriG.size() << endl;

	string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
	ifstream ain(AttFileName.c_str());

	if (!ain)
	{
		cout << "Fail to read " << AttFileName << "." << endl;
		return;
	}

	string aline;
	// unordered_map<int, set<int>> nodeAtt;

	while (getline(ain, aline))
	{
		string ver_s = aline.substr(0, aline.find("\t"));
		string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
		int ver = stof(ver_s); // string2float
		if (!oriG.count(ver))
			continue;
		int start = 0, end = 0;
		set<int> att;
		int w = 0;
		while ((end = att_s.find(",", start)) != string::npos)
		{
			w = stof(att_s.substr(start, end - start));
			att.insert(w);
			// auto it = att2vid.find(w);
			// if (it == att2vid.end())
			// {
			// 	att2vid.emplace(w, set<int>{ver});
			// }
			// else
			// {
			// 	att2vid[w].insert(ver);
			// }
			start = end + 1;
		}
		w = stof(att_s.substr(start));
		att.insert(w);
		// auto it = att2vid.find(w);
		// if (it == att2vid.end())
		// {
		// 	att2vid.emplace(w, set<int>{ver});
		// }
		// else
		// {
		// 	att2vid[w].insert(ver);
		// }

		id2att.emplace(ver, att);
	}
	ain.close();
	cout << "Loaded attribute successfully!" << endl;
	// findConnectedComponents(oriG, componentMap);
	string LocationFileName = "DataGraph/" + dataset + "/" + "user_location.txt";
	ifstream lin(LocationFileName.c_str());

	if (!lin)
	{
		cout << "No location information." << endl;
		// return;
	}
	else
	{
		string lline;
		while (getline(lin, lline))
		{
			if (lline.find('#') != string::npos)
				continue;
			size_t first_tab = lline.find("\t");
			string user_id = lline.substr(0, first_tab);
			size_t second_tab = lline.find("\t", first_tab + 1);
			string locs1 = lline.substr(first_tab + 1, second_tab - first_tab - 1);
			string locs2 = lline.substr(second_tab + 1);
			int uid = stof(user_id); // string2float
			if (!oriG.count(uid))
				continue;
			int loc1 = stof(locs1);
			int loc2 = stof(locs2);
			id2loc.emplace(uid, make_pair(loc1, loc2));
		}
	}

	lin.close();
	// remove on information vertices

	unordered_set<int> inter;
	if (id2loc.empty())
	{
		for (auto it : id2att)
		{
			inter.insert(it.first);
		}
	}
	else
	{
		for (auto it : id2att)
		{
			if (id2loc.count(it.first))
			{
				inter.insert(it.first);
			}
		}
	}

	std::vector<int> toRemove;
	for (const auto &it : oriG)
	{
		if (!inter.count(it.first))
		{
			toRemove.push_back(it.first);
		}
	}

	for (int id : toRemove)
	{
		
		for (int v : oriG[id])
		{
			oriG[v].erase(id);
		}
		if (id2att.count(id))
		{
			id2att.erase(id);
		}
		if (id2loc.count(id))
		{
			id2loc.erase(id);
		}
		oriG.erase(id);
	}
	for (auto it1 : id2att)
	{
		for (int w : it1.second)
		{
			auto it = att2vid.find(w);
			if (it == att2vid.end())
			{
				att2vid.emplace(w, set<int>{it1.first});
			}
			else
			{
				att2vid[w].insert(it1.first);
			}
		}
	}
	// unordered_set<int> visited;
	// unordered_set<int> remove;
	// for (auto it = oriG.begin(); it != oriG.end();)
	// {

	// 	if (it->second.size() < 2)
	// 	{
	// 		for (int v : it->second)
	// 		{
	// 			oriG[v].erase(it->first);
	// 			if (visited.count(v) && oriG[v].size() < 2)
	// 			{
	// 				remove.insert(v);
	// 			}
	// 		}
	// 		it = oriG.erase(it);
	// 	}
	// 	else
	// 	{
	// 		visited.insert(it->first);
	// 		++it;
	// 	}
	// }
	// queue<int> q;
	// for (int v : remove)
	// {
	// 	q.push(v);
	// }
	// while (!q.empty())
	// {
	// 	int v = q.front();
	// 	q.pop();

	// 	for (int u : oriG[v])
	// 	{
	// 		oriG[u].erase(v);
	// 		if (!remove.count(u) && oriG[u].size() < 2)
	// 		{
	// 			remove.insert(u);
	// 			q.push(u);
	// 		}
	// 	}
	// 	oriG.erase(v);
	// }
	cout << "Loaded dataset successfully!" << endl;
}
void KSubG(int k){
	unordered_set<int> visited;
	unordered_set<int> remove;
	for (auto it = oriG.begin(); it != oriG.end();)
	{

		if (it->second.size() < k-1)
		{
			for (int v : it->second)
			{
				oriG[v].erase(it->first);
				if (visited.count(v) && oriG[v].size() < k-1)
				{
					remove.insert(v);
				}
			}
			it = oriG.erase(it);
		}
		else
		{
			visited.insert(it->first);
			++it;
		}
	}
	queue<int> q;
	for (int v : remove)
	{
		q.push(v);
	}
	while (!q.empty())
	{
		int v = q.front();
		q.pop();

		for (int u : oriG[v])
		{
			oriG[u].erase(v);
			if (!remove.count(u) && oriG[u].size() < k-1)
			{
				remove.insert(u);
				q.push(u);
			}
		}
		oriG.erase(v);
	}
	for(int v : remove){
		if(id2loc.count(v)){
			id2loc.erase(v);
		}
		if(id2att.count(v)){
			id2att.erase(v);
		}
	}
}
void KMaxCS(set<int> &SCCResult, map<int, vector<Edge>, greater<int>> &rel2edge, unordered_map<int, set<int>> &rel2verts)
{
	DataGraph *relSubG = new DataGraph();
	auto rel = rel2edge.begin();
	while (rel != rel2edge.end())
	{
		for (auto &e : rel->second)
		{
			int src = e.first;
			int dst = e.second;
			relSubG->addEdgeNoMatinC(src, dst);
		}
		rel++;
	}
	TrussDecomposition *kt = new TrussDecomposition(relSubG);
	DataGraph *SubG = new DataGraph();
	for (auto &it : kt->k2edge[kt->kMax])
	{
		int src = relSubG->seq2id[it.first];
		int dst = relSubG->seq2id[it.second];
		SubG->addEdge(src, dst);
	}
	int minSize = INT16_MAX;
	for (auto &cc : SubG->CC)
	{
		if (cc.second.size() < minSize)
		{
			minSize = cc.first;
		}
	}
	unordered_map<int, unordered_map<int, int>> e2sups;
	SubG->support(e2sups, SubG->CC[minSize]);
	bool del = true;
	while (del)
	{
		del = false;
		for (int v : SubG->CC[minSize])
		{
			if (e2sups.count(v) && checkInd(v, kt->kMax, e2sups))
			{
				del = true;
			}
		}
	}
	// string s = "";
	// int eNum = 0;
	// for (auto &it : e2sups)
	// {
	// 	// s += to_string(SubG->seq2id[it.first]) + ",";
	// 	eNum += it.second.size();
	// }
	// eNum /= 2;
	// s += to_string(eNum);
	// cout<<"CC: "<<s<<endl;
	// SCCResult.push_back(s);
	for (auto &it : e2sups)
	{
		SCCResult.insert(SubG->seq2id[it.first]);
	}
}
void qv(set<int> &tempC, set<int> &qw, double &aaks, double &aapj, double &aapjk, double &aed, double &aapjkt, string type)
{
	if (tempC.empty())
		return;
	int vNum = tempC.size();

	double aks = 0;
	for (int v : tempC)
	{
		set<int> cap;
		set_intersection(id2att[v].begin(), id2att[v].end(), qw.begin(), qw.end(), inserter(cap, cap.begin()));
		int up = cap.size();
		int down = id2att[v].size() + qw.size() - up;
		aks += ((double)up / (double)down);
	}
	aks /= (double)vNum;

	double apj = 0;
	double apjk = 0;
	double apjkt = 0;
	vector<int> vv(tempC.begin(), tempC.end());
	for (int i = 0; i < vv.size(); i++)
	{
		for (int j = i + 1; j < vv.size(); j++)
		{
			int u = vv[i], v = vv[j];
			if (u == v)
				continue;

			if (type == "i")
			{
				set<int> inter;
				set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

				double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
				apj += in / un;
			}
			else
			{
				if (id2loc.count(u) && id2loc.count(v))
				{
					// apj += sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2));
					// apj += exp(-pow(id2loc[u].first - id2loc[v].first, 2) - pow(id2loc[u].second - id2loc[v].second, 2));
					apj += exp(-0.01 * (sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2))));
					cout << apj << endl;
				}
			}
			// set<int> cap;
			// set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(cap, cap.begin()));
			// int up = cap.size();
			// int down = id2att[u].size() + id2att[v].size() - up;
			// apj += ((double)up / (double)down);
			// set<int> capKw;
			// set_intersection(cap.begin(), cap.end(), qw.begin(), qw.end(), inserter(capKw, capKw.begin()));
			// int upKw = capKw.size();
			// apjkt += (double)upKw / (double)(cap.size() + qw.size() - upKw);
			// if (up == 0)
			// {
			// 	continue;
			// }
			// apjk += ((double)upKw / (double)up);
		}
	}
	apj /= (vNum * (vNum - 1) / 2);
	// apjk /= (vNum * (vNum - 1) / 2);
	// apjkt /= (vNum * (vNum - 1) / 2);

	double ed = 0;
	for (int u : tempC)
	{
		for (int v : oriG[u])
		{
			if (tempC.find(v) != tempC.end())
			{
				ed++;
			}
		}
	}
	ed /= (vNum * (vNum - 1));
	// s += to_string(eNum);
	// SCCResult.push_back(s);
	aaks += aks;
	aapj += apj;
	// aapjk += apjk;
	aed += ed;
	// aapjkt += apjkt;
}
void qv(set<int> &tempC, set<int> &qw, double &aaks, double &aapj, double &aapjk, double &aapjkt, string type)
{
	if (tempC.empty())
		return;
	int vNum = tempC.size();
	double aks = 0;
	for (int v : tempC)
	{
		set<int> cap;
		set_intersection(id2att[v].begin(), id2att[v].end(), qw.begin(), qw.end(), inserter(cap, cap.begin()));
		int up = cap.size();
		int down = id2att[v].size() + qw.size() - up;
		aks += ((double)up / (double)down);
	}
	aks /= (double)vNum;

	double apj = 0;
	double apjk = 0;
	double apjkt = 0;
	vector<int> vv(tempC.begin(), tempC.end());
	for (int i = 0; i < vv.size(); i++)
	{
		for (int j = i + 1; j < vv.size(); j++)
		{
			int u = vv[i], v = vv[j];
			if (u == v)
				continue;
			// set<int> cap;
			// set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(cap, cap.begin()));
			// int up = cap.size();

			// int down = id2att[u].size() + id2att[v].size() - up;
			// apj += ((double)up / (double)down);
			// set<int> capKw;
			// set_intersection(cap.begin(), cap.end(), qw.begin(), qw.end(), inserter(capKw, capKw.begin()));
			// int upKw = capKw.size();
			// apjkt += (double)upKw / (double)(cap.size() + qw.size() - upKw);
			// if (up == 0)
			// {
			// 	continue;
			// }
			if (type == "i")
			{
				set<int> inter;
				set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

				double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
				apj += in / un;
			}
			else
			{
				if (id2loc.count(u) && id2loc.count(v))
					// apj += sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2));
					// apj += exp(-pow(id2loc[u].first - id2loc[v].first, 2) - pow(id2loc[u].second - id2loc[v].second, 2));
					apj += exp(-0.01 * (sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2))));
			}
			// apjk += ((double)upKw / (double)up);
		}
	}
	apj /= (vNum * (vNum - 1) / 2);
	// apjk /= (vNum * (vNum - 1) / 2);
	// apjkt /= (vNum * (vNum - 1) / 2);

	aaks += aks;
	aapj += apj;
	// aapjk += apjk;
	// aapjkt += apjkt;
}
void readResult(string dataset, vector<set<int>> &result)
{
	string AttFileName = dataset + "_quality_result.txt";
	ifstream queryFile(AttFileName.c_str());
	string aline;
	double aed = 0;
	int notEmpty = 0;

	while (getline(queryFile, aline))
	{
		if (aline == "NULL")
		{
			result.push_back(set<int>());
			continue;
		}
		++notEmpty;
		int start = 0, end = 0;
		set<int> vertices;
		int w = 0;
		while ((end = aline.find(",", start)) != string::npos)
		{
			w = stof(aline.substr(start, end - start));
			vertices.insert(w);
			start = end + 1;
		}
		w = stof(aline.substr(start));
		vertices.insert(w);
		result.push_back(vertices);

		getline(queryFile, aline);
		double eNum = 0;
		std::istringstream iss(aline);
		std::string token;

		while (std::getline(iss, token, '\t'))
		{
			++eNum;
		}
		int vNum = vertices.size();
		eNum /= (vNum * (vNum - 1));
		aed += eNum;
	}
	queryFile.close();
	cout << "Average edge density: " << aed / notEmpty << endl;
}
// void LoadSimE(string &dataset)
// {
// 	string AttFileName = "DataGraph/" + dataset + "_simG.txt";
// 	ifstream ain(AttFileName.c_str());

// 	if (!ain)
// 	{
// 		cout << "Fail to read " << AttFileName << ", now recomputing." << endl;
// 	}
// 	else
// 	{
// 		int u, v;
// 		double w;
// 		while (ain >> u >> v >> w)
// 		{
// 			existV.insert(u);
// 			existV.insert(v);
// 			unordered_map<int, double> ew1;
// 			ew1[v] = w;
// 			auto ret1 = simE.emplace(u, ew1);
// 			if (!ret1.second)
// 			{
// 				simE[u][v] = w;
// 			}
// 			unordered_map<int, double> ew2;
// 			ew2[u] = w;
// 			auto ret2 = simE.emplace(v, ew2);
// 			if (!ret2.second)
// 			{
// 				simE[v][u] = w;
// 			}
// 		}
// 	}
// 	ain.close();
// }
// void Effectiveness(string &dataset, double &simThreshold)
// {
// 	vector<set<int>> kacResult;
// 	readResult("kac", kacResult);
// 	vector<set<int>> kcResult;
// 	readResult("kc", kcResult);
// 	string AttFileName = "DataGraph/" + dataset + "/" + "quality_query.txt";
// 	ifstream queryFile(AttFileName.c_str());

// 	int queryCnt;
// 	queryFile >> queryCnt;
// 	// vector<set<int>> queryAttrs;
// 	double aaks = 0;
// 	double aapj = 0;
// 	double aapjk = 0;
// 	double aapjkt = 0;
// 	double aed = 0;
// 	double aaksr = 0;
// 	double aapjr = 0;
// 	double aapjkr = 0;
// 	double aapjktr = 0;
// 	double aedr = 0;
// 	double aaks1 = 0;
// 	double aapj1 = 0;
// 	double aapjk1 = 0;
// 	double aapjkt1 = 0;
// 	double aed1 = 0;
// 	double aaks2 = 0;
// 	double aapj2 = 0;
// 	double aapjk2 = 0;
// 	double aapjkt2 = 0;
// 	double aed2 = 0;
// 	double aaks3 = 0;
// 	double aapj3 = 0;
// 	double aapjk3 = 0;
// 	double aapjkt3 = 0;
// 	double aed3 = 0;

// 	int myNum = 0;
// 	int sccNum = 0;
// 	int rkacNum = 0;
// 	int kcNum = 0;
// 	for (int i = 0; i < queryCnt; i++)
// 	{
// 		string queryAttrStr;
// 		set<int> qw;
// 		queryFile >> queryAttrStr;
// 		int pos1 = 0, pos2 = 0;
// 		while ((pos2 = queryAttrStr.find(",", pos1)) != string::npos)
// 		{
// 			qw.insert(stof(queryAttrStr.substr(pos1, pos2 - pos1)));
// 			pos1 = pos2 + 1;
// 		}
// 		qw.insert(stof(queryAttrStr.substr(pos1)));
// 		// queryAttrs.push_back(queryAttr);
// 		map<int, vector<Edge>, greater<int>> rel2edge;
// 		unordered_map<int, set<int>> rel2verts;
// 		getKWSubG(dataset, qw, rel2edge, rel2verts);
// 		set<int> SCCResult;
// 		KMaxCS(SCCResult, rel2edge, rel2verts);

// 		// int conQw = rel2edge.begin()->first;
// 		// if (conQw > 1)
// 		// 	cout << "qw number: " << conQw << i << endl;
// 		auto start = std::chrono::high_resolution_clock::now();
// 		map<int, unordered_map<int, double>> simE;
// 		int relmax = getSimEs(simThreshold, dataset, rel2verts, simE, i);
// 		// for (auto it = rel2edge.begin(); it != rel2edge.end();)
// 		// {
// 		// 	if (it->first < relmax)
// 		// 	{
// 		// 		it = rel2edge.erase(it);
// 		// 	}
// 		// 	else
// 		// 	{
// 		// 		++it;
// 		// 	}
// 		// }

// 		auto end = std::chrono::high_resolution_clock::now();
// 		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
// 		std::cout << "getSimEs : " << duration.count() << " milliseconds." << std::endl;
// 		vector<SkyGroupCand> result;

// 		start = std::chrono::high_resolution_clock::now();
// 		Best(rel2edge, result, simE, rel2verts);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
// 		std::cout << "Best: " << duration.count() << " milliseconds." << std::endl;

// 		set<int> tempC1;
// 		// set<int> tempC2;
// 		int rmax = 0;
// 		int kmax = 0;
// 		// double rhomax = 0;
// 		for (auto &r : result)
// 		{
// 			if (r.rel >= rmax)
// 			{
// 				rmax = r.rel;
// 				if (r.k > kmax)
// 				{
// 					tempC1 = r.vertices;
// 					kmax = r.k;
// 				}
// 			}
// 			// if (r.rho > rhomax)
// 			// {
// 			// 	tempC2 = r.vertices;
// 			// 	rhomax = r.rho;
// 			// }
// 		}
// 		if (!result.empty())
// 			++myNum;
// 		// string s = "";
// 		// for (int v : tempC)
// 		// {
// 		// 	s += to_string(v) + ",";
// 		// }
// 		qv(tempC1, qw, aaks, aapj, aapjk, aed, aapjkt);
// 		// qv(tempC2, qw, aaksr, aapjr, aapjkr, aedr, aapjktr);

// 		if (!SCCResult.empty())
// 			++sccNum;
// 		qv(SCCResult, qw, aaks1, aapj1, aapjk1, aed1, aapjkt1);
// 		if (!kacResult[i].empty())
// 			++rkacNum;
// 		qv(kacResult[i], qw, aaks2, aapj2, aapjk2, aapjkt2);
// 		if (!kcResult[i].empty())
// 			++kcNum;
// 		qv(kcResult[i], qw, aaks3, aapj3, aapjk3, aapjkt3);
// 	}
// 	queryFile.close();
// 	aaks /= myNum;
// 	aapj /= myNum;
// 	aapjk /= myNum;
// 	aed /= myNum;
// 	aapjkt /= myNum;
// 	cout << "Average aks: " << aaks << endl;
// 	cout << "Average apj: " << aapj << endl;
// 	cout << "Average apjk: " << aapjk << endl;
// 	cout << "Average apjkt: " << aapjkt << endl;
// 	cout << "Average ed: " << aed << endl;
// 	cout << "Average sr: " << myNum / 100.0 << endl;

// 	// aaksr /= myNum;
// 	// aapjr /= myNum;
// 	// aapjkr /= myNum;
// 	// aedr /= myNum;
// 	// aapjktr /= myNum;
// 	// cout << "Average aksr: " << aaksr<< endl;
// 	// cout << "Average apjr: " << aapjr << endl;
// 	// cout << "Average apjkr: " << aapjkr << endl;
// 	// cout << "Average apjktr: " << aapjktr << endl;
// 	// cout << "Average edr: " << aedr << endl;

// 	aaks1 /= sccNum;
// 	aapj1 /= sccNum;
// 	aapjk1 /= sccNum;
// 	aed1 /= sccNum;
// 	aapjkt1 /= sccNum;
// 	cout << "Average aks1: " << aaks1 << endl;
// 	cout << "Average apj1: " << aapj1 << endl;
// 	cout << "Average apjk1: " << aapjk1 << endl;
// 	cout << "Average apjkt1: " << aapjkt1 << endl;
// 	cout << "Average ed1: " << aed1 << endl;
// 	cout << "Average sccsr: " << sccNum / 100.0 << endl;

// 	aaks2 /= rkacNum;
// 	aapj2 /= rkacNum;
// 	aapjk2 /= rkacNum;
// 	aapjkt2 /= rkacNum;
// 	cout << "Average aks2: " << aaks2 << endl;
// 	cout << "Average apj2: " << aapj2 << endl;
// 	cout << "Average apjk2: " << aapjk2 << endl;
// 	cout << "Average apjkt2: " << aapjkt2 << endl;
// 	cout << "Average kacsr: " << rkacNum / 100.0 << endl;

// 	aaks3 /= kcNum;
// 	aapj3 /= kcNum;
// 	aapjk3 /= kcNum;
// 	aapjkt3 /= kcNum;
// 	cout << "Average aks3: " << aaks3 << endl;
// 	cout << "Average apj3: " << aapj3 << endl;
// 	cout << "Average apjk3: " << aapjk3 << endl;
// 	cout << "Average apjkt3: " << aapjkt3 << endl;
// 	cout << "Average kcsr: " << kcNum / 100.0 << endl;
// }
void Effectiveness(string &dataset, double &simThreshold, string type)
{
	vector<set<int>> kacResult;
	readResult("kac", kacResult);
	vector<set<int>> kcResult;
	readResult("kc", kcResult);
	string AttFileName = "DataGraph/" + dataset + "/" + "quality_query.txt";
	ifstream queryFile(AttFileName.c_str());

	int queryCnt;
	queryFile >> queryCnt;
	// vector<set<int>> queryAttrs;
	double aaks = 0;
	double aapj = 0;
	double aapjk = 0;
	double aapjkt = 0;
	double aed = 0;
	double aaksr = 0;
	double aapjr = 0;
	double aapjkr = 0;
	double aapjktr = 0;
	double aedr = 0;
	double aaks1 = 0;
	double aapj1 = 0;
	double aapjk1 = 0;
	double aapjkt1 = 0;
	double aed1 = 0;
	double aaks2 = 0;
	double aapj2 = 0;
	double aapjk2 = 0;
	double aapjkt2 = 0;
	double aed2 = 0;
	double aaks3 = 0;
	double aapj3 = 0;
	double aapjk3 = 0;
	double aapjkt3 = 0;
	double aed3 = 0;

	int myNum = 0;
	int sccNum = 0;
	int rkacNum = 0;
	int kcNum = 0;
	for (int i = 0; i < queryCnt; i++)
	{
		string queryAttrStr;
		set<int> qw;
		queryFile >> queryAttrStr;
		int pos1 = 0, pos2 = 0;
		while ((pos2 = queryAttrStr.find(",", pos1)) != string::npos)
		{
			qw.insert(stof(queryAttrStr.substr(pos1, pos2 - pos1)));
			pos1 = pos2 + 1;
		}
		qw.insert(stof(queryAttrStr.substr(pos1)));
		// queryAttrs.push_back(queryAttr);
		map<int, vector<Edge>, greater<int>> rel2edge;
		unordered_map<int, set<int>> rel2verts;
		getKWSubG(dataset, qw, rel2edge, rel2verts);
		set<int> SCCResult;
		KMaxCS(SCCResult, rel2edge, rel2verts);

		int relmax = rel2edge.begin()->first;

		map<int, unordered_map<int, double>> simE;
		getMaxRelSimEs(relmax, dataset, rel2verts, simE, type);
		for (auto it = rel2edge.begin(); it != rel2edge.end();)
		{
			if (it->first < relmax)
			{
				it = rel2edge.erase(it);
			}
			else
			{
				++it;
			}
		}

		vector<SkyGroupCand> result;
		if (i == 21)
		{
			bool flag = false;
		}
		Best(rel2edge, result, simE, rel2verts);

		set<int> tempC1;
		// set<int> tempC2;
		int rmax = 0;
		int kmax = 0;
		// double rhomax = 0;
		for (auto &r : result)
		{
			if (r.rel >= rmax)
			{
				rmax = r.rel;
				if (r.k > kmax)
				{
					tempC1 = r.vertices;
					kmax = r.k;
				}
			}
			// if (r.rho > rhomax)
			// {
			// 	tempC2 = r.vertices;
			// 	rhomax = r.rho;
			// }
		}
		if (!result.empty())
			++myNum;
		// string s = "";
		// for (int v : tempC)
		// {
		// 	s += to_string(v) + ",";
		// }
		qv(tempC1, qw, aaks, aapj, aapjk, aed, aapjkt, type);
		// qv(tempC2, qw, aaksr, aapjr, aapjkr, aedr, aapjktr);

		if (!SCCResult.empty())
			++sccNum;
		qv(SCCResult, qw, aaks1, aapj1, aapjk1, aed1, aapjkt1, type);
		if (!kacResult[i].empty())
			++rkacNum;
		qv(kacResult[i], qw, aaks2, aapj2, aapjk2, aapjkt2, type);
		if (!kcResult[i].empty())
			++kcNum;
		qv(kcResult[i], qw, aaks3, aapj3, aapjk3, aapjkt3, type);
	}
	queryFile.close();
	aaks /= myNum;
	aapj /= myNum;
	aed /= myNum;
	cout << "Average aks: " << aaks << endl;
	cout << "Average apj: " << aapj << endl;
	cout << "Average ed: " << aed << endl;
	cout << "Average sr: " << myNum / 100.0 << endl;

	aaks1 /= sccNum;
	aapj1 /= sccNum;
	aed1 /= sccNum;
	cout << "Average aks1: " << aaks1 << endl;
	cout << "Average apj1: " << aapj1 << endl;
	cout << "Average ed1: " << aed1 << endl;
	cout << "Average sccsr: " << sccNum / 100.0 << endl;

	aaks2 /= rkacNum;
	aapj2 /= rkacNum;
	cout << "Average aks2: " << aaks2 << endl;
	cout << "Average apj2: " << aapj2 << endl;
	cout << "Average kacsr: " << rkacNum / 100.0 << endl;

	aaks3 /= kcNum;
	aapj3 /= kcNum;
	cout << "Average aks3: " << aaks3 << endl;
	cout << "Average apj3: " << aapj3 << endl;
	cout << "Average kcsr: " << kcNum / 100.0 << endl;
}

void Approximate(string dataset, double simThreshold, string type)
{
	string AttFileName = "DataGraph/" + dataset + "/" + "quality_query.txt";
	ifstream queryFile(AttFileName.c_str());

	double basicA = 0;
	double improveA = 0;
	double bestA = 0;
	int basicNum = 0;
	int improveNum = 0;
	int bestNum = 0;
	int queryCnt;
	double ba_im = 0.0;
	double ba_be = 0.0;
	double im_be = 0.0;
	queryFile >> queryCnt;
	for (int i = 0; i < queryCnt; i++)
	{
		cout << "query: " << i << endl;
		string queryAttrStr;
		set<int> qw;
		queryFile >> queryAttrStr;
		int pos1 = 0, pos2 = 0;
		while ((pos2 = queryAttrStr.find(",", pos1)) != string::npos)
		{
			qw.insert(stof(queryAttrStr.substr(pos1, pos2 - pos1)));
			pos1 = pos2 + 1;
		}
		qw.insert(stof(queryAttrStr.substr(pos1)));
		// queryAttrs.push_back(queryAttr);
		map<int, vector<Edge>, greater<int>> rel2edge;
		unordered_map<int, set<int>> rel2verts;
		getKWSubG(dataset, qw, rel2edge, rel2verts);
		int relmax = rel2edge.begin()->first;

		map<int, unordered_map<int, double>> simE;
		getMaxRelSimEs(relmax, dataset, rel2verts, simE, type);

		// getSimEs(simThreshold, dataset, rel2verts, simE, i);
		// 找到rel最大且可构成图
		int64_t n = simE.size();
		vector<int> ei;
		vector<int> ej;
		vector<double> ew;
		unordered_map<int, int> visited;
		unordered_map<int, int> seq2id;
		int id = 0;

		for (auto &uv : simE)
		{
			int u = uv.first;
			for (auto &vw : uv.second)
			{
				int v = vw.first;
				double w = vw.second;
				if (!visited.count(u))
				{
					visited.emplace(u, id);
					seq2id.emplace(id, u);
					id++;
				}
				int u20 = visited[u];
				if (!visited.count(v))
				{
					visited.emplace(v, id);
					seq2id.emplace(id, v);
					id++;
				}
				int v20 = visited[v];
				ei.push_back(u20);
				ej.push_back(v20);
				ew.push_back(w);
			}
		}
		int64_t m = ei.size();
		if (m == 0)
			continue;
		for (auto it = rel2edge.begin(); it != rel2edge.end();)
		{
			if (it->first < relmax)
			{
				it = rel2edge.erase(it);
			}
			else
			{
				++it;
			}
		}
		vector<int> output;

		double mop = densest_subgraph(n, m, ei, ej, ew, output) / 2;
		if (i == 41)
		{
			bool flag = false;
		}
		cout << "mop: " << mop << endl;
		int minsup = INT16_MAX;
		for (int i = 0; i < output.size(); i++)
		{
			output[i] = seq2id[output[i]];
		}
		for (int u : output)
		{

			for (int v : oriG[u])
			{
				if (v > u && find(output.begin(), output.end(), v) != output.end())
				{
					int sup = 0;
					for (int w : oriG[v])
					{
						if (find(output.begin(), output.end(), w) != output.end() && oriG[u].count(w))
						{
							sup++;
						}
					}
					minsup = min(minsup, sup);
				}
			}
		}
		// find minsup+2-truss densest subgraph
		vector<SkyGroupCand> result;
		// if(i == 8){
		// 	bool flag = false;
		// }
		Basic1(minsup + 1, rel2edge, result, simE, rel2verts);
		double basic = 0.0;
		if (!result.empty())
		{
			basic = result[0].rho;
			basicA += mop / basic;
			result.clear();
			basicNum++;
			cout << "BasicRho: " << basic << endl;
			cout << "BasicK: " << result[0].k << endl;
			cout << "BasicRel: " << result[0].rel << endl;
		}
		Improve1(minsup + 1, rel2edge, result, simE, rel2verts);
		double improve = 0.0;
		if (!result.empty())
		{
			improve = result[0].rho;
			improveA += mop / improve;
			ba_im += basic / improve;
			result.clear();
			improveNum++;
			cout << "ImproveRho: " << improve << endl;
			cout << "ImproveK: " << result[0].k << endl;
			cout << "ImproveRel: " << result[0].rel << endl;
		}
		Best1(minsup + 1, rel2edge, result, simE, rel2verts);
		if (!result.empty())
		{
			double best = result[0].rho;
			bestA += mop / best;
			ba_be += basic / best;
			im_be += improve / best;
			result.clear();
			bestNum++;
			cout << "BestRho: " << best << endl;
			cout << "BestK: " << result[0].k << endl;
			cout << "BestRel: " << result[0].rel << endl;
		}
	}
	cout << "Average basic: " << basicA / basicNum << endl;
	cout << "Average improve: " << improveA / improveNum << endl;
	cout << "Average best: " << bestA / bestNum << endl;
	cout << "Average ba_im: " << ba_im / basicNum << endl;
	cout << "Average ba_be: " << ba_be / basicNum << endl;
	cout << "Average im_be: " << im_be / improveNum << endl;
	queryFile.close();
}

void Efficiency(int kUp, string dataset, double simThreshold, string type)
{
	string AttFileName = "DataGraph/" + dataset + "/" + "quality_query.txt";
	ifstream queryFile(AttFileName.c_str());

	int queryCnt;
	queryFile >> queryCnt;
	vector<set<int>> qws;
	for (int i = 0; i < queryCnt; i++)
	{
		string queryAttrStr;
		set<int> qw;
		queryFile >> queryAttrStr;
		int pos1 = 0, pos2 = 0;
		while ((pos2 = queryAttrStr.find(",", pos1)) != string::npos)
		{
			qw.insert(stof(queryAttrStr.substr(pos1, pos2 - pos1)));
			pos1 = pos2 + 1;
		}
		qw.insert(stof(queryAttrStr.substr(pos1)));
		qws.push_back(qw);
	}
	queryFile.close();
	// int k = 3;
	cout << "simThreshold: " << simThreshold << endl;
	for (int k = 3; k <= kUp; k++)
	{
		cout << "k: " << k << endl;
		double basicAvgTime = 0.0;
		double improveAvgTime = 0.0;
		double bestAvgTime = 0.0;
		double readQwAvgTime = 0.0;
		double readSimAvgTime = 0.0;
		for (int i = 0; i < queryCnt; i++)
		{
			cout << "query: " << i << endl;
			set<int> &qw = qws[i];
			// queryAttrs.push_back(queryAttr);
			map<int, vector<Edge>, greater<int>> rel2edge;
			unordered_map<int, set<int>> rel2verts;
			auto startr = std::chrono::high_resolution_clock::now();
			getKWSubG(dataset, qw, rel2edge, rel2verts);
			auto endr = std::chrono::high_resolution_clock::now();
			auto durationr = std::chrono::duration_cast<std::chrono::milliseconds>(endr - startr);
			readQwAvgTime += durationr.count();
			map<int, unordered_map<int, double>> simE;
			auto starts = std::chrono::high_resolution_clock::now();
			getSimEs(dataset, rel2verts, simE, type, simThreshold);
			auto ends = std::chrono::high_resolution_clock::now();
			auto durations = std::chrono::duration_cast<std::chrono::milliseconds>(ends - starts);
			readSimAvgTime += durations.count();

			// find minsup+2-truss densest subgraph
			vector<SkyGroupCand> result;

			auto start1 = std::chrono::high_resolution_clock::now();
			Basic(k, rel2edge, result, simE, rel2verts);
			auto end1 = std::chrono::high_resolution_clock::now();
			auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
			basicAvgTime += duration1.count();
			// if (!result.empty())
			// 	successNum++;
			result.clear();
			auto start2 = std::chrono::high_resolution_clock::now();
			Improve(k, rel2edge, result, simE, rel2verts);
			auto end2 = std::chrono::high_resolution_clock::now();
			auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
			improveAvgTime += duration2.count();
			result.clear();
			auto start3 = std::chrono::high_resolution_clock::now();
			Best(k, rel2edge, result, simE, rel2verts);
			auto end3 = std::chrono::high_resolution_clock::now();
			auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3);
			bestAvgTime += duration3.count();
			result.clear();
		}
		cout << "Average basic: " << basicAvgTime / queryCnt << " milliseconds" << endl;
		cout << "Average improve: " << improveAvgTime / queryCnt << " milliseconds" << endl;
		cout << "Average best: " << bestAvgTime / queryCnt << " milliseconds" << endl;
		cout << "Average basic+qw: " << (basicAvgTime + readQwAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average improve+qw: " << (improveAvgTime + readQwAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average best+qw: " << (bestAvgTime + readQwAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average basic+qw+sim: " << (basicAvgTime + readQwAvgTime + readSimAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average improve+qw+sim: " << (improveAvgTime + readQwAvgTime + readSimAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average best+qw+sim: " << (bestAvgTime + readQwAvgTime + readSimAvgTime) / queryCnt << " milliseconds" << endl;
	}
}

void QEfficiency(int k, string dataset, double simThreshold, string type)
{
	string AttFileName = "DataGraph/" + dataset + "/" + "quality_query.txt";
	ifstream queryFile(AttFileName.c_str());

	int queryCnt;
	queryFile >> queryCnt;
	vector<vector<int>> qws;
	for (int i = 0; i < queryCnt; i++)
	{
		string queryAttrStr;
		vector<int> qw;
		queryFile >> queryAttrStr;
		int pos1 = 0, pos2 = 0;
		while ((pos2 = queryAttrStr.find(",", pos1)) != string::npos)
		{
			qw.push_back(stof(queryAttrStr.substr(pos1, pos2 - pos1)));
			pos1 = pos2 + 1;
		}
		qw.push_back(stof(queryAttrStr.substr(pos1)));
		qws.push_back(qw);
	}
	queryFile.close();
	// int k = 3;
	KSubG(k);
	for (int j = 0; j < 5; j++)
	{

		cout << "qw: " << j + 1 << endl;
		double basicAvgTime = 0.0;
		double improveAvgTime = 0.0;
		double bestAvgTime = 0.0;
		double readQwAvgTime = 0.0;
		double readSimAvgTime = 0.0;
		for (int i = 0; i < queryCnt; i++)
		{
			cout << "query: " << i << endl;
			vector<int> &qwv = qws[i];
			cout << "qw: ";
			set<int> qw;
			for (int a = 0; a <= j; a++)
			{
				qw.insert(qwv[a]);
				cout << qwv[a] << " ";
			}
			cout << endl;
			// queryAttrs.push_back(queryAttr);
			map<int, vector<Edge>, greater<int>> rel2edge;
			unordered_map<int, set<int>> rel2verts;
			auto startr = std::chrono::high_resolution_clock::now();
			getKWSubG(dataset, qw, rel2edge, rel2verts);
			auto endr = std::chrono::high_resolution_clock::now();
			auto durationr = std::chrono::duration_cast<std::chrono::milliseconds>(endr - startr);
			readQwAvgTime += durationr.count();
			// getKWSubG(dataset, qw, rel2edge, rel2verts);

			map<int, unordered_map<int, double>> simE;
			// getSimEs(dataset, rel2verts, simE, type, simThreshold);
			auto starts = std::chrono::high_resolution_clock::now();
			getSimEs(dataset, rel2verts, simE, type, simThreshold);
			auto ends = std::chrono::high_resolution_clock::now();
			auto durations = std::chrono::duration_cast<std::chrono::milliseconds>(ends - starts);
			readSimAvgTime += durations.count();

			// find minsup+2-truss densest subgraph
			vector<SkyGroupCand> result;

			auto start1 = std::chrono::high_resolution_clock::now();
			Basic(k, rel2edge, result, simE, rel2verts);
			auto end1 = std::chrono::high_resolution_clock::now();
			auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
			basicAvgTime += duration1.count();
			// if (!result.empty())
			// 	successNum++;
			result.clear();
			auto start2 = std::chrono::high_resolution_clock::now();
			Improve(k, rel2edge, result, simE, rel2verts);
			auto end2 = std::chrono::high_resolution_clock::now();
			auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
			improveAvgTime += duration2.count();
			result.clear();
			auto start3 = std::chrono::high_resolution_clock::now();
			Best(k, rel2edge, result, simE, rel2verts);
			auto end3 = std::chrono::high_resolution_clock::now();
			auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3);
			bestAvgTime += duration3.count();
			result.clear();
		}
		cout << "Average basic: " << basicAvgTime / queryCnt << " milliseconds" << endl;
		cout << "Average improve: " << improveAvgTime / queryCnt << " milliseconds" << endl;
		cout << "Average best: " << bestAvgTime / queryCnt << " milliseconds" << endl;
		cout << "Average basic+qw: " << (basicAvgTime + readQwAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average improve+qw: " << (improveAvgTime + readQwAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average best+qw: " << (bestAvgTime + readQwAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average basic+qw+sim: " << (basicAvgTime + readQwAvgTime + readSimAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average improve+qw+sim: " << (improveAvgTime + readQwAvgTime + readSimAvgTime) / queryCnt << " milliseconds" << endl;
		cout << "Average best+qw+sim: " << (bestAvgTime + readQwAvgTime + readSimAvgTime) / queryCnt << " milliseconds" << endl;
	}
}
// 计算所有边的相似度，统计相似度的分布
// void CalESim(string type, string dataset)
// {
// 	unordered_map<double, int> similarityCount;
// 	for (const auto &[u, neighbors] : oriG)
// 	{
// 		for (int v : neighbors)
// 		{
// 			if (u < v)
// 			{ // 确保每条边只计算一次
// 				double w = 0;
// 				if (type == "i")
// 				{
// 					if (id2att.count(u) && id2att.count(v))
// 					{
// 						set<int> inter;
// 						set_intersection(id2att[u].begin(), id2att[u].end(), id2att[v].begin(), id2att[v].end(), inserter(inter, inter.begin()));

// 						double in = inter.size(), un = id2att[u].size() + id2att[v].size() - inter.size();
// 						w = in / un;
// 					}
// 				}
// 				else
// 				{
// 					if (id2loc.count(u) && id2loc.count(v))
// 						// w = sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2));
// 						// w = exp(-pow(id2loc[u].first - id2loc[v].first, 2) - pow(id2loc[u].second - id2loc[v].second, 2));
// 						w = exp(-0.01 * (sqrt(pow(id2loc[u].first - id2loc[v].first, 2) + pow(id2loc[u].second - id2loc[v].second, 2))));
// 				}

// 				similarityCount[w]++;
// 			}
// 		}
// 	}

// 	// 输出每个相似度及其对应的边的数量到文件
// 	string outputFileName = "DataGraph/" + dataset + "/" + "sim2enum.txt";
// 	ofstream outFile(outputFileName);

// 	for (const auto &[similarity, count] : similarityCount)
// 	{
// 		outFile << similarity << "\t" << count << endl;
// 	}

// 	outFile.close();
// }
int main(int argc, char *argv[])
{
	string dataset = argv[1];
	double simThreshold = stod(argv[2]);
	string method = argv[3];
	string dataType = argv[4];
	string ks = argv[5];
	int k = stoi(ks);
	// LoadData(k, dataset);
	LoadData(dataset);
	// CalESim(dataType, dataset);
	if (method == "0")
		Effectiveness(dataset, simThreshold, dataType);
	else if (method == "1")
		Approximate(dataset, simThreshold, dataType);
	else if (method == "2")
	{
		Efficiency(k, dataset, simThreshold, dataType);
	}
	else if(method == "3")
	{
		QEfficiency(k, dataset, simThreshold, dataType);
	}
	else{
		vector<double> simThresholds = {0.2, 0.4, 0.6, 0.8};
		for(int i = 0; i < simThresholds.size(); i++){
			Efficiency(k, dataset, simThresholds[i], dataType);
		}
	}
	// string AttFileName = "DataGraph/" + dataset + "/" + "quality_query.txt";
	// ifstream queryFile(AttFileName.c_str());

	// int queryCnt;
	// queryFile >> queryCnt;
	// for (int i = 0; i < queryCnt; i++)
	// {
	// 	string queryAttrStr;
	// 	set<int> qw;
	// 	queryFile >> queryAttrStr;
	// 	int pos1 = 0, pos2 = 0;
	// 	while ((pos2 = queryAttrStr.find(",", pos1)) != string::npos)
	// 	{
	// 		qw.insert(stof(queryAttrStr.substr(pos1, pos2 - pos1)));
	// 		pos1 = pos2 + 1;
	// 	}
	// 	qw.insert(stof(queryAttrStr.substr(pos1)));
	// 	// queryAttrs.push_back(queryAttr);
	// 	map<int, vector<Edge>, greater<int>> rel2edge;
	// 	unordered_map<int, set<int>> rel2verts;
	// 	getKWSubG(dataset, qw, rel2edge, rel2verts);

	// 	// int conQw = rel2edge.begin()->first;
	// 	// if (conQw > 1)
	// 	// 	cout << "qw number: " << conQw << i << endl;
	// 	auto start = std::chrono::high_resolution_clock::now();
	// 	map<int, unordered_map<int, double>> simE;
	// 	int relmax = getSimEs(simThreshold, dataset, rel2verts, simE, i);
	// 	// for (auto it = rel2edge.begin(); it != rel2edge.end();)
	// 	// {
	// 	// 	if (it->first < relmax)
	// 	// 	{
	// 	// 		it = rel2edge.erase(it);
	// 	// 	}
	// 	// 	else
	// 	// 	{
	// 	// 		++it;
	// 	// 	}
	// 	// }

	// 	auto end = std::chrono::high_resolution_clock::now();
	// 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	// 	std::cout << "getSimEs : " << duration.count() << " milliseconds." << std::endl;
	// 	vector<SkyGroupCand> result;

	// 	start = std::chrono::high_resolution_clock::now();
	// 	if (method == "0")
	// 		Basic(rel2edge, result, simE, rel2verts);
	// 	else if (method == "1")
	// 		Improve(rel2edge, result, simE, rel2verts);
	// 	else
	// 		Best(rel2edge, result, simE, rel2verts);
	// 	end = std::chrono::high_resolution_clock::now();
	// 	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	// 	std::cout << method << " : " << duration.count() << " milliseconds." << std::endl;
	// }
	// queryFile.close();
	return 0;
}