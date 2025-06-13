#include "DataGraph.h"
#include <algorithm>
#include <memory>

// bool compVertex(int i, int j)
// {
// 	return AdjList[i].size() < AdjList[j].size() || (AdjList[i].size() == AdjList[j].size() && i < j);
// 	// return 1;
// }

// void makeSet(shared_ptr<UNode> &x)
// {
// 	x->parent = x;
// 	x->rank = 0;
// }

// shared_ptr<UNode> find(shared_ptr<UNode> &x)
// {
// 	if (x->parent != x)
// 	{
// 		x->parent = find(x->parent); // make the parent of x point to the root
// 	}
// 	return x->parent;
// }
void dfs(int v, vector<vector<int>> &adj, vector<bool> &visited)
{
	visited[v] = true;
	for (int u : adj[v])
	{
		if (!visited[u])
		{
			dfs(u, adj, visited);
		}
	}
}
bool isConnected(vector<vector<int>> &adj)
{
	vector<bool> visited(adj.size(), false); // 记录顶点是否访问过

	// 从顶点 0 开始 DFS
	dfs(0, adj, visited);

	// 检查所有顶点是否都被访问过
	for (bool v : visited)
	{
		if (!v)
		{
			return false; // 如果有顶点没被访问过，图不连通
		}
	}
	return true;
}
void DataGraph::unite(shared_ptr<UNode> &x, shared_ptr<UNode> &y) // union
{
	shared_ptr<UNode> xRoot = find(x);
	shared_ptr<UNode> yRoot = find(y);

	// Compare the underlying values rather than the pointers
	if (xRoot->value == yRoot->value)
	{
		return;
	}

	// x and y are not already in the same set. Merge them.
	if (xRoot->rank < yRoot->rank)
	{
		xRoot->parent = yRoot;
		for (int node : CC[xRoot->value])
		{
			// CC[yRoot->value].push_back(node);
			CC[yRoot->value].insert(node);
		}
		CC.erase(xRoot->value);
	}
	else if (xRoot->rank > yRoot->rank)
	{
		yRoot->parent = xRoot;
		for (int node : CC[yRoot->value])
		{
			// CC[xRoot->value].push_back(node);
			CC[xRoot->value].insert(node);
		}
		CC.erase(yRoot->value);
	}
	else
	{
		yRoot->parent = xRoot;
		xRoot->rank = xRoot->rank + 1;
		for (int node : CC[yRoot->value])
		{
			// CC[xRoot->value].push_back(node);
			CC[xRoot->value].insert(node);
		}
		CC.erase(yRoot->value);
	}
}
DataGraph::DataGraph() : vNum(0), eNum(0) {}
// DataGraph::DataGraph(string &fileName, int &N, int &M)
// {
// 	// vector<vector<int>>().swap(AdjList);
// 	// unordered_map<int, int>().swap(id2seq);
// 	// unordered_map<int, int>().swap(seq2id);
// 	// seq2att.clear();

// 	///////////////////////////////////////////////////////////////////////////////////////////
// 	string StructFileName = "DataGraph/" + fileName + ".txt";
// 	ifstream sin(StructFileName.c_str());

// 	if (!sin)
// 	{
// 		cout << "Fail to read " << StructFileName << "." << endl;
// 		return;
// 	}

// 	N = 0, M = 0;

// 	string sline;

// 	while (getline(sin, sline))
// 	{
// 		if (sline.find('#') != string::npos)
// 			continue;
// 		string src_s = sline.substr(0, sline.find("\t"));
// 		string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
// 		int src = stof(src_s); // string2float
// 		int dst = stof(dst_s);
// 		if (src != dst)
// 		{
// 			// cout<<src<<" "<<dst<<endl;
// 			unordered_map<int, int>::iterator it; // unordered_map is the hashtable
// 			int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

// 			it = id2seq.find(src);
// 			if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
// 			{
// 				seq1 = AdjList.size();
// 				id2seq[src] = seq1;
// 				seq2id[seq1] = src;
// 				vector<int> newVertex;
// 				AdjList.push_back(newVertex);
// 				// seq2att[seq1] = nodeAtt[src];
// 			}
// 			else
// 			{
// 				seq1 = it->second;
// 			}

// 			it = id2seq.find(dst);
// 			if (it == id2seq.end())
// 			{
// 				seq2 = AdjList.size();
// 				id2seq[dst] = seq2;
// 				seq2id[seq2] = dst;
// 				vector<int> newVertex;
// 				AdjList.push_back(newVertex);
// 				// seq2att[seq2] = nodeAtt[dst];
// 			}
// 			else
// 			{
// 				seq2 = it->second;
// 			}

// 			int flag;

// 			flag = 0;
// 			for (int i = 0; i < AdjList[seq1].size(); i++)
// 			{
// 				if (AdjList[seq1][i] == seq2)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

// 			flag = 0;
// 			for (int i = 0; i < AdjList[seq2].size(); i++)
// 			{
// 				if (AdjList[seq2][i] == seq1)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				AdjList[seq2].push_back(seq1), M++;
// 			// cout<<M<<endl;
// 		}
// 	}
// 	sin.close();
// 	N = id2seq.size();
// 	cout << "vertices = " << N << " ; edges = " << M << endl;
// 	// write file
// 	//  for (int i = 0; i < N; ++i)
// 	//  	verSort.push_back(i);
// 	//  sort(verSort.begin(), verSort.end(), [this](int i, int j)
// 	//  	 { return AdjList[i].size() < AdjList[j].size() || (AdjList[i].size() == AdjList[j].size() && i < j); });
// 	////////////////////////////////////////////////////////////////////////////////
// 	string AttFileName = "DataGraph/" + fileName + "_attribute.txt";
// 	ifstream ain(AttFileName.c_str());

// 	if (!ain)
// 	{
// 		cout << "Fail to read " << AttFileName << "." << endl;
// 		return;
// 	}

// 	string aline;
// 	unordered_map<int, set<int>> nodeAtt;

// 	while (getline(ain, aline))
// 	{
// 		string ver_s = aline.substr(0, aline.find("\t"));
// 		string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
// 		int ver = stof(ver_s); // string2float

// 		int start = 0, end = 0;
// 		set<int> att;
// 		while ((end = att_s.find(",", start)) != string::npos)
// 		{
// 			att.insert(stof(att_s.substr(start, end - start)));
// 			start = end + 1;
// 		}
// 		att.insert(stof(att_s.substr(start)));
// 		int seq = id2seq[ver];
// 		if (seq2att.find(seq) == seq2att.end())
// 		{
// 			seq2att[seq] = att;
// 		}
// 		// if (nodeAtt.find(ver) == nodeAtt.end())
// 		// {
// 		// 	nodeAtt[ver] = att;
// 		// }
// 	}
// 	ain.close();
// 	// write file
// 	///////////////////////////////////////////////////////////////////////////////////////

// 	// vector<int> keys;
// 	// for (const auto &pair : seq2att)
// 	// {
// 	// 	unordered_map<int, double> adj;
// 	// 	SimAdjList.push_back(adj);
// 	// 	keys.push_back(pair.first);
// 	// }
// 	// sort(keys.begin(), keys.end());

// 	// for (size_t i = 0; i < keys.size(); ++i)
// 	// {
// 	// 	for (size_t j = i + 1; j < keys.size(); ++j)
// 	// 	{
// 	// 		int key1 = keys[i];
// 	// 		int key2 = keys[j];

// 	// 		const set<int> &set1 = seq2att.at(key1);
// 	// 		const set<int> &set2 = seq2att.at(key2);

// 	// 		vector<int> intersection;
// 	// 		set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
// 	// 						 back_inserter(intersection));

// 	// 		double intersection_count = intersection.size();

// 	// 		if (intersection_count > 0)
// 	// 		{
// 	// 			SimAdjList[key1].insert(make_pair(key2, intersection_count));

// 	// 			SimAdjList[key2].insert(make_pair(key1, intersection_count));
// 	// 		}
// 	// 	}
// 	// }
// }
DataGraph::DataGraph(string &fileName, int &N, int &M)
{
	// vector<vector<int>>().swap(AdjList);
	// unordered_map<int, int>().swap(id2seq);
	// unordered_map<int, int>().swap(seq2id);
	// seq2att.clear();

	///////////////////////////////////////////////////////////////////////////////////////////
	string StructFileName = "DataGraph/" + fileName + ".txt";
	ifstream sin(StructFileName.c_str());

	if (!sin)
	{
		cout << "Fail to read " << StructFileName << "." << endl;
		return;
	}

	N = 0, M = 0;

	string sline;

	while (getline(sin, sline))
	{
		if (sline.find('#') != string::npos)
			continue;
		string src_s = sline.substr(0, sline.find("\t"));
		string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
		int src = stof(src_s); // string2float
		int dst = stof(dst_s);
		if (src != dst)
		{
			// cout<<src<<" "<<dst<<endl;
			unordered_map<int, int>::iterator it; // unordered_map is the hashtable
			int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

			it = id2seq.find(src);
			if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
			{
				seq1 = AdjList.size();
				id2seq[src] = seq1;
				seq2id[seq1] = src;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq1] = nodeAtt[src];
			}
			else
			{
				seq1 = it->second;
			}

			it = id2seq.find(dst);
			if (it == id2seq.end())
			{
				seq2 = AdjList.size();
				id2seq[dst] = seq2;
				seq2id[seq2] = dst;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq2] = nodeAtt[dst];
			}
			else
			{
				seq2 = it->second;
			}

			int flag;

			flag = 0;
			for (int i = 0; i < AdjList[seq1].size(); i++)
			{
				if (AdjList[seq1][i] == seq2)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

			flag = 0;
			for (int i = 0; i < AdjList[seq2].size(); i++)
			{
				if (AdjList[seq2][i] == seq1)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				AdjList[seq2].push_back(seq1), M++;
			// cout<<M<<endl;
		}
	}
	sin.close();
	N = id2seq.size();
	cout << "vertices = " << N << " ; edges = " << M << endl;
	// write file
	//  for (int i = 0; i < N; ++i)
	//  	verSort.push_back(i);
	//  sort(verSort.begin(), verSort.end(), [this](int i, int j)
	//  	 { return AdjList[i].size() < AdjList[j].size() || (AdjList[i].size() == AdjList[j].size() && i < j); });
	////////////////////////////////////////////////////////////////////////////////
	string AttFileName = "DataGraph/" + fileName + "_attribute.txt";
	ifstream ain(AttFileName.c_str());

	if (!ain)
	{
		cout << "Fail to read " << AttFileName << "." << endl;
		return;
	}

	string aline;
	unordered_map<int, set<int>> nodeAtt;

	while (getline(ain, aline))
	{
		string ver_s = aline.substr(0, aline.find("\t"));
		string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
		int ver = stof(ver_s); // string2float

		int start = 0, end = 0;
		set<int> att;
		while ((end = att_s.find(",", start)) != string::npos)
		{
			att.insert(stof(att_s.substr(start, end - start)));
			start = end + 1;
		}
		att.insert(stof(att_s.substr(start)));
		if (id2att.find(ver) == id2att.end())
		{
			id2att[ver] = att;
		}
		// if (nodeAtt.find(ver) == nodeAtt.end())
		// {
		// 	nodeAtt[ver] = att;
		// }
	}
	ain.close();
}
DataGraph::DataGraph(vector<pair<int, int>> &edges)
{
	// vector<vector<int>>().swap(AdjList);
	// unordered_map<int, int>().swap(id2seq);
	// unordered_map<int, int>().swap(seq2id);
	// seq2att.clear();

	///////////////////////////////////////////////////////////////////////////////////////////

	for (Edge e : edges)
	{
		int src = e.first;
		int dst = e.second;
		if (src != dst)
		{
			// cout<<src<<" "<<dst<<endl;
			unordered_map<int, int>::iterator it; // unordered_map is the hashtable
			int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

			it = id2seq.find(src);
			if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
			{
				seq1 = AdjList.size();
				id2seq[src] = seq1;
				seq2id[seq1] = src;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq1] = nodeAtt[src];
			}
			else
			{
				seq1 = it->second;
			}

			it = id2seq.find(dst);
			if (it == id2seq.end())
			{
				seq2 = AdjList.size();
				id2seq[dst] = seq2;
				seq2id[seq2] = dst;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq2] = nodeAtt[dst];
			}
			else
			{
				seq2 = it->second;
			}

			int flag;

			flag = 0;
			for (int i = 0; i < AdjList[seq1].size(); i++)
			{
				if (AdjList[seq1][i] == seq2)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

			flag = 0;
			for (int i = 0; i < AdjList[seq2].size(); i++)
			{
				if (AdjList[seq2][i] == seq1)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				AdjList[seq2].push_back(seq1);
			// cout<<M<<endl;
		}
	}
}
//

// DataGraph::DataGraph(unordered_map<int, int> &vert2rel, DataGraph *dg)
// {
// 	// vector<vector<int>>().swap(AdjList);
// 	// unordered_map<int, int>().swap(id2seq);
// 	// unordered_map<int, int>().swap(seq2id);
// 	// seq2att.clear();

// 	///////////////////////////////////////////////////////////////////////////////////////////
// 	vector<int> verts;
// 	for (auto &vr : vert2rel)
// 	{
// 		verts.push_back(vr.first);
// 	}
// 	sort(verts.begin(), verts.end());
// 	for (auto &v : verts)
// 	{
// 		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
// 		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst
// 		int id = dg->seq2id[v];
// 		it = id2seq.find(id);
// 		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
// 		{
// 			seq1 = AdjList.size();
// 			id2seq[id] = seq1;
// 			seq2id[seq1] = id;
// 			id2att[id] = dg->id2att[id];
// 			vector<int> newVertex;
// 			AdjList.push_back(newVertex);
// 			int rel = vert2rel[v];
// 			seq2rel[seq1] = rel;

// 			if (rel2seqs.find(rel) == rel2seqs.end())
// 			{
// 				set<int> seqs = {seq1};
// 				rel2seqs[rel] = seqs;
// 			}
// 			else
// 			{
// 				rel2seqs[rel].insert(seq1);
// 			}
// 		}
// 		else
// 		{
// 			seq1 = it->second;
// 		}

// 		vector<int> nei = dg->AdjList[v];
// 		sort(nei.begin(), nei.end());
// 		vector<int> com;
// 		set_intersection(nei.begin(), nei.end(), verts.begin(), verts.end(), back_inserter(com));
// 		for (auto &n : com)
// 		{
// 			id = dg->seq2id[n];
// 			it = id2seq.find(id);
// 			if (it == id2seq.end())
// 			{
// 				seq2 = AdjList.size();
// 				id2seq[id] = seq2;
// 				seq2id[seq2] = id;
// 				id2att[id] = dg->id2att[id];
// 				vector<int> newVertex;
// 				AdjList.push_back(newVertex);
// 				int rel = vert2rel[n];
// 				seq2rel[seq2] = rel;

// 				if (rel2seqs.find(rel) == rel2seqs.end())
// 				{
// 					set<int> seqs = {seq2};
// 					rel2seqs[rel] = seqs;
// 				}
// 				else
// 				{
// 					rel2seqs[rel].insert(seq2);
// 				}
// 			}
// 			else
// 			{
// 				seq2 = it->second;
// 			}

// 			int flag;

// 			flag = 0;
// 			for (int i = 0; i < AdjList[seq1].size(); i++)
// 			{
// 				if (AdjList[seq1][i] == seq2)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

// 			flag = 0;
// 			for (int i = 0; i < AdjList[seq2].size(); i++)
// 			{
// 				if (AdjList[seq2][i] == seq1)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				AdjList[seq2].push_back(seq1);
// 		}
// 	}
// }
// DataGraph *DataGraph::getSubG(int rel)
// {
// 	DataGraph *subG = new DataGraph();

// 	set<int> verts;
// 	map<int, set<int>>::reverse_iterator it;
// 	for (it = this->rel2seqs.rbegin(); it != this->rel2seqs.rend(); it++)
// 	{
// 		if (it->first < rel)
// 			break;
// 		verts.insert(it->second.begin(), it->second.end());
// 	}
// 	for (auto &v : verts)
// 	{
// 		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
// 		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst
// 		int id = this->seq2id[v];
// 		it = subG->id2seq.find(id);
// 		if (it == subG->id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
// 		{
// 			seq1 = subG->AdjList.size();
// 			subG->id2seq[id] = seq1;
// 			subG->seq2id[seq1] = id;
// 			vector<int> newVertex;
// 			subG->AdjList.push_back(newVertex);
// 		}
// 		else
// 		{
// 			seq1 = it->second;
// 		}

// 		vector<int> nei = this->AdjList[v];
// 		sort(nei.begin(), nei.end());
// 		vector<int> com;
// 		set_intersection(nei.begin(), nei.end(), verts.begin(), verts.end(), back_inserter(com));
// 		for (auto &n : com)
// 		{
// 			id = this->seq2id[n];
// 			it = subG->id2seq.find(id);
// 			if (it == subG->id2seq.end())
// 			{
// 				seq2 = subG->AdjList.size();
// 				subG->id2seq[id] = seq2;
// 				subG->seq2id[seq2] = id;
// 				vector<int> newVertex;
// 				subG->AdjList.push_back(newVertex);
// 			}
// 			else
// 			{
// 				seq2 = it->second;
// 			}

// 			int flag;

// 			flag = 0;
// 			for (int i = 0; i < subG->AdjList[seq1].size(); i++)
// 			{
// 				if (subG->AdjList[seq1][i] == seq2)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				subG->AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

// 			flag = 0;
// 			for (int i = 0; i < subG->AdjList[seq2].size(); i++)
// 			{
// 				if (subG->AdjList[seq2][i] == seq1)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				subG->AdjList[seq2].push_back(seq1);
// 		}
// 	}
// 	return subG;
// }

// DataGraph *DataGraph::getSubG(set<pair<int, int>> &edges)
// {

// 	DataGraph *subG = new DataGraph();

// 	// N = 0, M = 0;
// 	int n = 0;

// 	string sline;
// 	map<int, shared_ptr<UNode>> idUFMap;

// 	for (auto e : edges)
// 	{
// 		int src = e.first; // string2float
// 		int dst = e.second;
// 		if (src != dst)
// 		{
// 			// cout<<src<<" "<<dst<<endl;
// 			unordered_map<int, int>::iterator it;	// unordered_map is the hashtable
// 			int seq1, seq2, id = this->seq2id[src]; // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst
// 			it = subG->id2seq.find(id);
// 			if (it == subG->id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
// 			{
// 				n++;
// 				seq1 = subG->AdjList.size();
// 				shared_ptr<UNode> unode = make_shared<UNode>(seq1);
// 				makeSet(unode);
// 				idUFMap[seq1] = unode;
// 				subG->id2seq[id] = seq1;
// 				subG->seq2id[seq1] = id;
// 				vector<int> newVertex;
// 				subG->AdjList.push_back(newVertex);
// 			}
// 			else
// 			{
// 				seq1 = it->second;
// 			}
// 			id = this->seq2id[dst];
// 			it = subG->id2seq.find(id);
// 			if (it == subG->id2seq.end())
// 			{
// 				n++;
// 				seq2 = subG->AdjList.size();
// 				shared_ptr<UNode> unode = make_shared<UNode>(seq2);
// 				makeSet(unode);
// 				idUFMap[seq2] = unode;
// 				subG->id2seq[id] = seq2;
// 				subG->seq2id[seq2] = id;
// 				vector<int> newVertex;
// 				subG->AdjList.push_back(newVertex);
// 			}
// 			else
// 			{
// 				seq2 = it->second;
// 			}
// 			unite(idUFMap[seq1], idUFMap[seq2], this->CC);

// 			int flag;

// 			flag = 0;
// 			for (int i = 0; i < subG->AdjList[seq1].size(); i++)
// 			{
// 				if (subG->AdjList[seq1][i] == seq2)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				subG->AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

// 			flag = 0;
// 			for (int i = 0; i < subG->AdjList[seq2].size(); i++)
// 			{
// 				if (subG->AdjList[seq2][i] == seq1)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				subG->AdjList[seq2].push_back(seq1);
// 			// M++;
// 			// cout<<M<<endl;
// 		}
// 	}

// 	return subG;
// }
void DataGraph::getSimGs(double simThreshold, map<int, unordered_map<int, double>> &simG, unordered_map<int, set<int>> &id2att, int ccId)
{

	unordered_map<int, set<int>> kw2vert;
	for (auto &c : CC[ccId])
	{
		for (int kw : id2att[seq2id[c]])
		{
			kw2vert[kw].insert(c);
		}
	}

	// vector<int> verts = cc.second;
	// sort(verts.begin(), verts.end());
	set<int> verts = CC[ccId];
	for (auto &c : verts)
	{

		set<int> nei;
		for (int kw : id2att[seq2id[c]])
		{
			nei.insert(kw2vert[kw].begin(), kw2vert[kw].end());
		}

		for (int v : nei)
		{
			if (v <= c)
				continue;
			set<int> inter;
			set_intersection(id2att[seq2id[c]].begin(), id2att[seq2id[c]].end(), id2att[seq2id[v]].begin(), id2att[seq2id[v]].end(), inserter(inter, inter.begin()));
			set<int> unin;
			set_union(id2att[seq2id[c]].begin(), id2att[seq2id[c]].end(), id2att[seq2id[v]].begin(), id2att[seq2id[v]].end(), inserter(unin, unin.begin()));
			// Weight w(inter.size(), unin.size());
			// Weight thd(0,5);
			double in = inter.size(), un = unin.size();
			double w = in / un;
			if (inter.size() > 0 && w > simThreshold)
			{
				unordered_map<int, double> ew1;
				ew1[v] = w;
				auto ret1 = simG.emplace(c, ew1);
				if (!ret1.second)
				{
					simG[c][v] = w;
				}
				unordered_map<int, double> ew2;
				ew2[c] = w;
				auto ret2 = simG.emplace(v, ew2);
				if (!ret2.second)
				{
					simG[v][c] = w;
				}
			}
		}
	}
}
// unordered_map<int, map<int, unordered_map<int, double>>> &DataGraph::getSimGs(unordered_map<int, set<int>> &id2att)
// {
// 	unordered_map<int, unordered_map<int, set<int>>> kw2verts;

// 	for (auto &cc : CC)
// 	{
// 		unordered_map<int, set<int>> kw2vert;
// 		for (auto &c : cc.second)
// 		{
// 			for (int kw : id2att[seq2id[c]])
// 			{
// 				kw2vert[kw].insert(c);
// 			}
// 		}
// 		kw2verts[cc.first] = kw2vert;
// 	}
// 	for (auto &cc : CC)
// 	{

// 		// vector<int> verts = cc.second;
// 		// sort(verts.begin(), verts.end());
// 		set<int> verts = cc.second;
// 		map<int, unordered_map<int, double>> g;
// 		for (auto &c : verts)
// 		{

// 			set<int> nei;
// 			for (int kw : id2att[seq2id[c]])
// 			{
// 				nei.insert(kw2verts[cc.first][kw].begin(), kw2verts[cc.first][kw].end());
// 			}

// 			for (int v : nei)
// 			{
// 				if (v <= c)
// 					continue;
// 				set<int> inter;
// 				set_intersection(id2att[seq2id[c]].begin(), id2att[seq2id[c]].end(), id2att[seq2id[v]].begin(), id2att[seq2id[v]].end(), inserter(inter, inter.begin()));
// 				set<int> unin;
// 				set_union(id2att[seq2id[c]].begin(), id2att[seq2id[c]].end(), id2att[seq2id[v]].begin(), id2att[seq2id[v]].end(), inserter(unin, unin.begin()));
// 				// Weight w(inter.size(), unin.size());
// 				// Weight thd(0,5);
// 				double in = inter.size(), un = unin.size();
// 				double w = in / un;
// 				if (inter.size() > 0 && w > simThreshold)
// 				{
// 					unordered_map<int, double> ew1;
// 					ew1[v] = w;
// 					auto ret1 = g.emplace(c, ew1);
// 					if (!ret1.second)
// 					{
// 						g[c][v] = w;
// 					}
// 					unordered_map<int, double> ew2;
// 					ew2[c] = w;
// 					auto ret2 = g.emplace(v, ew2);
// 					if (!ret2.second)
// 					{
// 						g[v][c] = w;
// 					}
// 				}
// 			}
// 		}
// 		SimG[cc.first] = g;
// 	}
// 	return SimG;
// }
void DataGraph::addEdges(vector<Edge> &newEs, vector<pair<int, int>> &edges)
{
	for (auto &e : edges)
	{
		int src = e.first; // string2float
		int dst = e.second;
		if (src != dst)
		{
			// cout<<src<<" "<<dst<<endl;
			unordered_map<int, int>::iterator it; // unordered_map is the hashtable
			int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

			it = id2seq.find(src);
			if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
			{
				seq1 = AdjList.size();
				id2seq[src] = seq1;
				seq2id[seq1] = src;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq1] = nodeAtt[src];
			}
			else
			{
				seq1 = it->second;
			}

			it = id2seq.find(dst);
			if (it == id2seq.end())
			{
				seq2 = AdjList.size();
				id2seq[dst] = seq2;
				seq2id[seq2] = dst;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq2] = nodeAtt[dst];
			}
			else
			{
				seq2 = it->second;
			}

			int flag;

			flag = 0;
			for (int i = 0; i < AdjList[seq1].size(); i++)
			{
				if (AdjList[seq1][i] == seq2)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
			{
				AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded
				newEs.push_back(Edge(min(seq1, seq2), max(seq1, seq2)));
			}

			flag = 0;
			for (int i = 0; i < AdjList[seq2].size(); i++)
			{
				if (AdjList[seq2][i] == seq1)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				AdjList[seq2].push_back(seq1);
		}
	}
}
void DataGraph::addEdge(int src, int dst)
{
	if (src != dst)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

		it = id2seq.find(src);
		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			seq1 = AdjList.size();
			shared_ptr<UNode> unode = make_shared<UNode>(seq1);
			makeSet(unode);
			idUFMap.emplace(seq1, unode);
			// CC[seq1].push_back(seq1);
			CC[seq1].insert(seq1);
			id2seq[src] = seq1;
			seq2id[seq1] = src;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq1] = nodeAtt[src];
		}
		else
		{
			seq1 = it->second;
		}

		it = id2seq.find(dst);
		if (it == id2seq.end())
		{
			seq2 = AdjList.size();
			shared_ptr<UNode> unode = make_shared<UNode>(seq2);
			makeSet(unode);
			idUFMap.emplace(seq2, unode);
			// CC[seq2].push_back(seq2);
			CC[seq2].insert(seq2);
			id2seq[dst] = seq2;
			seq2id[seq2] = dst;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq2] = nodeAtt[dst];
		}
		else
		{
			seq2 = it->second;
		}
		unite(idUFMap[seq1], idUFMap[seq2]);
		int flag;

		flag = 0;
		for (int i = 0; i < AdjList[seq1].size(); i++)
		{
			if (AdjList[seq1][i] == seq2)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

		flag = 0;
		for (int i = 0; i < AdjList[seq2].size(); i++)
		{
			if (AdjList[seq2][i] == seq1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq2].push_back(seq1);
	}
}
void DataGraph::addEdgeNoMatinC(int src, int dst)
{
	if (src != dst)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

		it = id2seq.find(src);
		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			seq1 = AdjList.size();
			id2seq[src] = seq1;
			seq2id[seq1] = src;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq1] = nodeAtt[src];
		}
		else
		{
			seq1 = it->second;
		}

		it = id2seq.find(dst);
		if (it == id2seq.end())
		{
			seq2 = AdjList.size();
			id2seq[dst] = seq2;
			seq2id[seq2] = dst;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq2] = nodeAtt[dst];
		}
		else
		{
			seq2 = it->second;
		}
		int flag;

		flag = 0;
		for (int i = 0; i < AdjList[seq1].size(); i++)
		{
			if (AdjList[seq1][i] == seq2)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

		flag = 0;
		for (int i = 0; i < AdjList[seq2].size(); i++)
		{
			if (AdjList[seq2][i] == seq1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq2].push_back(seq1);
	}
}
void DataGraph::addEdgesAndMaintainC(vector<pair<int, int>> &edges)
{
	for (auto &e : edges)
	{
		int src = e.first; // string2float
		int dst = e.second;
		if (src != dst)
		{
			// cout<<src<<" "<<dst<<endl;
			unordered_map<int, int>::iterator it; // unordered_map is the hashtable
			int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

			it = id2seq.find(src);
			if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
			{
				seq1 = AdjList.size();
				shared_ptr<UNode> unode = make_shared<UNode>(seq1);
				makeSet(unode);
				idUFMap.emplace(seq1, unode);
				CC[seq1].insert(seq1);
				id2seq[src] = seq1;
				seq2id[seq1] = src;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq1] = nodeAtt[src];
			}
			else
			{
				seq1 = it->second;
			}

			it = id2seq.find(dst);
			if (it == id2seq.end())
			{
				seq2 = AdjList.size();
				shared_ptr<UNode> unode = make_shared<UNode>(seq2);
				makeSet(unode);
				idUFMap.emplace(seq2, unode);
				CC[seq2].insert(seq2);
				id2seq[dst] = seq2;
				seq2id[seq2] = dst;
				vector<int> newVertex;
				AdjList.push_back(newVertex);
				// seq2att[seq2] = nodeAtt[dst];
			}
			else
			{
				seq2 = it->second;
			}
			unite(idUFMap[seq1], idUFMap[seq2]);
			int flag;

			flag = 0;
			for (int i = 0; i < AdjList[seq1].size(); i++)
			{
				if (AdjList[seq1][i] == seq2)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
			{
				AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded
			}

			flag = 0;
			for (int i = 0; i < AdjList[seq2].size(); i++)
			{
				if (AdjList[seq2][i] == seq1)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				AdjList[seq2].push_back(seq1);
		}
	}
}

void DataGraph::addEdgeAndMainConnect(int src, int dst, map<int, unordered_map<int, double>> &simE, unordered_map<int, unordered_map<int, double>> &simG)
{
	if (src != dst)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

		it = id2seq.find(src);
		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			seq1 = AdjList.size();
			shared_ptr<UNode> unode = make_shared<UNode>(seq1);
			makeSet(unode);
			idUFMap.emplace(seq1, unode);
			// CC[seq1].push_back(seq1);
			CC[seq1].insert(seq1);
			id2seq[src] = seq1;
			seq2id[seq1] = src;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq1] = nodeAtt[src];

			simG.emplace(seq1, unordered_map<int, double>());
			for (auto &t : simE[src])
			{
				auto tp = id2seq.find(t.first);
				if (tp != id2seq.end())
				{
					double w = simE[src][t.first];
					simG[seq1][tp->second] = w;
					simG[tp->second][seq1] = w;
				}
			}
		}
		else
		{
			seq1 = it->second;
		}

		it = id2seq.find(dst);
		if (it == id2seq.end())
		{
			seq2 = AdjList.size();
			shared_ptr<UNode> unode = make_shared<UNode>(seq2);
			makeSet(unode);
			idUFMap.emplace(seq2, unode);
			// CC[seq2].push_back(seq2);
			CC[seq2].insert(seq2);
			id2seq[dst] = seq2;
			seq2id[seq2] = dst;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq2] = nodeAtt[dst];

			simG.emplace(seq2, unordered_map<int, double>());
			for (auto &t : simE[dst])
			{
				auto tp = id2seq.find(t.first);
				if (tp != id2seq.end())
				{
					double w = simE[dst][t.first];
					simG[seq2][tp->second] = w;
					simG[tp->second][seq2] = w;
				}
			}
		}
		else
		{
			seq2 = it->second;
		}
		unite(idUFMap[seq1], idUFMap[seq2]);
		int flag;

		flag = 0;
		for (int i = 0; i < AdjList[seq1].size(); i++)
		{
			if (AdjList[seq1][i] == seq2)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

		flag = 0;
		for (int i = 0; i < AdjList[seq2].size(); i++)
		{
			if (AdjList[seq2][i] == seq1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq2].push_back(seq1);
	}
}
void DataGraph::addEdgeAndMainConnect(vector<Edge> &newEs, int src, int dst, map<int, unordered_map<int, double>> &simE, unordered_map<int, unordered_map<int, double>> &simG)
{
	if (src != dst)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

		it = id2seq.find(src);
		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			seq1 = AdjList.size();
			shared_ptr<UNode> unode = make_shared<UNode>(seq1);
			makeSet(unode);
			idUFMap.emplace(seq1, unode);
			// CC[seq1].push_back(seq1);
			CC[seq1].insert(seq1);
			id2seq[src] = seq1;
			seq2id[seq1] = src;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq1] = nodeAtt[src];

			simG.emplace(seq1, unordered_map<int, double>());
			for (auto &t : simE[src])
			{
				auto tp = id2seq.find(t.first);
				if (tp != id2seq.end())
				{
					double w = simE[src][t.first];
					simG[seq1][tp->second] = w;
					simG[tp->second][seq1] = w;
				}
			}
		}
		else
		{
			seq1 = it->second;
		}

		it = id2seq.find(dst);
		if (it == id2seq.end())
		{
			seq2 = AdjList.size();
			shared_ptr<UNode> unode = make_shared<UNode>(seq2);
			makeSet(unode);
			idUFMap.emplace(seq2, unode);
			// CC[seq2].push_back(seq2);
			CC[seq2].insert(seq2);
			id2seq[dst] = seq2;
			seq2id[seq2] = dst;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq2] = nodeAtt[dst];

			simG.emplace(seq2, unordered_map<int, double>());
			for (auto &t : simE[dst])
			{
				auto tp = id2seq.find(t.first);
				if (tp != id2seq.end())
				{
					double w = simE[dst][t.first];
					simG[seq2][tp->second] = w;
					simG[tp->second][seq2] = w;
				}
			}
		}
		else
		{
			seq2 = it->second;
		}
		unite(idUFMap[seq1], idUFMap[seq2]);
		int flag;

		flag = 0;
		for (int i = 0; i < AdjList[seq1].size(); i++)
		{
			if (AdjList[seq1][i] == seq2)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded
			newEs.push_back(Edge(min(seq1, seq2), max(seq1, seq2)));
		}

		flag = 0;
		for (int i = 0; i < AdjList[seq2].size(); i++)
		{
			if (AdjList[seq2][i] == seq1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq2].push_back(seq1);
	}
}
void DataGraph::addEdge(Edge &e)
{

	int src = e.first; // string2float
	int dst = e.second;
	if (src != dst)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

		it = id2seq.find(src);
		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			seq1 = AdjList.size();
			id2seq[src] = seq1;
			seq2id[seq1] = src;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq1] = nodeAtt[src];
		}
		else
		{
			seq1 = it->second;
		}

		it = id2seq.find(dst);
		if (it == id2seq.end())
		{
			seq2 = AdjList.size();
			id2seq[dst] = seq2;
			seq2id[seq2] = dst;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq2] = nodeAtt[dst];
		}
		else
		{
			seq2 = it->second;
		}

		int flag;

		flag = 0;
		for (int i = 0; i < AdjList[seq1].size(); i++)
		{
			if (AdjList[seq1][i] == seq2)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

		flag = 0;
		for (int i = 0; i < AdjList[seq2].size(); i++)
		{
			if (AdjList[seq2][i] == seq1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq2].push_back(seq1);
	}
}
DataGraph *DataGraph::getSubG(vector<int> &vertices)
{

	DataGraph *subG = new DataGraph();

	// N = 0, M = 0;

	string sline;
	// map<int, shared_ptr<UNode>> idUFMap;
	sort(vertices.begin(), vertices.end());
	for (auto v : vertices)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2, id = this->seq2id[v]; // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst
		it = subG->id2seq.find(id);
		if (it == subG->id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			subG->vNum++;
			seq1 = subG->AdjList.size();
			// shared_ptr<UNode> unode = make_shared<UNode>(seq1);
			// makeSet(unode);
			// idUFMap[seq1] = unode;
			subG->id2seq[id] = seq1;
			subG->seq2id[seq1] = id;
			vector<int> newVertex;
			subG->AdjList.push_back(newVertex);
		}
		else
		{
			seq1 = it->second;
		}
		vector<int> nei;
		sort(this->AdjList[v].begin(), this->AdjList[v].end());
		set_intersection(this->AdjList[v].begin(), this->AdjList[v].end(), vertices.begin(), vertices.end(), back_inserter(nei));
		for (int u : nei)
		{
			id = this->seq2id[u];
			it = subG->id2seq.find(id);
			if (it == subG->id2seq.end())
			{
				subG->vNum++;

				seq2 = subG->AdjList.size();
				// shared_ptr<UNode> unode = make_shared<UNode>(seq2);
				// makeSet(unode);
				// idUFMap[seq2] = unode;
				subG->id2seq[id] = seq2;
				subG->seq2id[seq2] = id;
				vector<int> newVertex;
				subG->AdjList.push_back(newVertex);
			}
			else
			{
				seq2 = it->second;
			}
			// unite(idUFMap[seq1], idUFMap[seq2]);
			if (seq1 > seq2)
				continue;

			int flag;

			flag = 0;
			for (int i = 0; i < subG->AdjList[seq1].size(); i++)
			{
				if (subG->AdjList[seq1][i] == seq2)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
			{
				subG->AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded
				subG->eNum++;
			}

			flag = 0;
			for (int i = 0; i < subG->AdjList[seq2].size(); i++)
			{
				if (subG->AdjList[seq2][i] == seq1)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				subG->AdjList[seq2].push_back(seq1);
		}

		// M++;
		// cout<<M<<endl;
	}
	// unordered_map<int, vector<int>> components;
	// for (int i = 0; i < vertices.size(); ++i)
	// {
	// 	int root = find(idUFMap[i])->value;
	// 	subG->CC[root].push_back(i);
	// }

	return subG;
}
void DataGraph::getCC(unordered_map<int, set<int>> &cc, set<int> &vertices)
{
	map<int, shared_ptr<UNode>> idUFMap;
	// unordered_set<int> visit;
	for (auto v : vertices)
	{
		shared_ptr<UNode> unode = make_shared<UNode>(v);
		makeSet(unode);
		idUFMap[v] = unode;
	}
	for (auto v : vertices)
	{
		vector<int> nei;
		sort(this->AdjList[v].begin(), this->AdjList[v].end());
		set_intersection(this->AdjList[v].begin(), this->AdjList[v].end(), vertices.begin(), vertices.end(), back_inserter(nei));
		for (int u : nei)
		{
			if (v > u)
				continue;
			unite(idUFMap[v], idUFMap[u]);
		}
	}
	for (int i : vertices)
	{
		int root = find(idUFMap[i])->value;
		auto it = cc.find(root);
		if (it == cc.end())
			cc.emplace(root, set<int>{i});
		else
			cc[root].insert(i);
	}
}
// void DataGraph::getSubG(unordered_map<int, set<int>> &cc, unordered_map<int, unordered_map<int, int>> &edge2sup, set<int> &vertices)
// {
// 	map<int, shared_ptr<UNode>> idUFMap;
// 	for (auto v : vertices)
// 	{
// 		auto it = edge2sup.find(v);
// 		if (it == edge2sup.end())
// 		{
// 			shared_ptr<UNode> unode = make_shared<UNode>(v);
// 			makeSet(unode);
// 			idUFMap[v] = unode;
// 			edge2sup.emplace(v, unordered_map<int, int>());
// 		}
// 		vector<int> nei;
// 		sort(this->AdjList[v].begin(), this->AdjList[v].end());
// 		set_intersection(this->AdjList[v].begin(), this->AdjList[v].end(), vertices.begin(), vertices.end(), back_inserter(nei));
// 		for (int u : nei)
// 		{
// 			if (v > u)
// 				continue;
// 			auto it2 = edge2sup.find(u);
// 			if (it2 == edge2sup.end())
// 			{
// 				shared_ptr<UNode> unode = make_shared<UNode>(u);
// 				makeSet(unode);
// 				idUFMap[u] = unode;

// 				edge2sup.emplace(v, unordered_map<int, int>());
// 			}
// 			unite(idUFMap[v], idUFMap[u]);
// 			edge2sup[u][v] = -1;
// 			edge2sup[v][u] = -1;
// 		}
// 		for (int i : vertices)
// 		{
// 			int root = find(idUFMap[i])->value;
// 			auto it = cc.find(root);
// 			if (it == cc.end())
// 				cc.emplace(root, set<int>{i});
// 			else
// 				cc[root].insert(i);
// 		}
// 	}
// }
// void DataGraph::getSubG(DataGraph *subG, set<int> &vertices)
// {

// 	// N = 0, M = 0;

// 	string sline;
// 	map<int, shared_ptr<UNode>> idUFMap;
// 	for (auto v : vertices)
// 	{
// 		// cout<<src<<" "<<dst<<endl;
// 		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
// 		int seq1, seq2, id = this->seq2id[v]; // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst
// 		it = subG->id2seq.find(id);
// 		if (it == subG->id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
// 		{
// 			subG->vNum++;
// 			seq1 = subG->AdjList.size();
// 			shared_ptr<UNode> unode = make_shared<UNode>(seq1);
// 			makeSet(unode);
// 			idUFMap[seq1] = unode;
// 			subG->id2seq[id] = seq1;
// 			subG->seq2id[seq1] = id;
// 			vector<int> newVertex;
// 			subG->AdjList.push_back(newVertex);
// 		}
// 		else
// 		{
// 			seq1 = it->second;
// 		}
// 		vector<int> nei;
// 		sort(this->AdjList[v].begin(), this->AdjList[v].end());
// 		set_intersection(this->AdjList[v].begin(), this->AdjList[v].end(), vertices.begin(), vertices.end(), back_inserter(nei));
// 		for (int u : nei)
// 		{
// 			id = this->seq2id[u];
// 			it = subG->id2seq.find(id);
// 			if (it == subG->id2seq.end())
// 			{
// 				subG->vNum++;

// 				seq2 = subG->AdjList.size();
// 				shared_ptr<UNode> unode = make_shared<UNode>(seq2);
// 				makeSet(unode);
// 				idUFMap[seq2] = unode;
// 				subG->id2seq[id] = seq2;
// 				subG->seq2id[seq2] = id;
// 				vector<int> newVertex;
// 				subG->AdjList.push_back(newVertex);
// 			}
// 			else
// 			{
// 				seq2 = it->second;
// 			}
// 			unite(idUFMap[seq1], idUFMap[seq2]);
// 			if (seq1 > seq2)
// 				continue;

// 			int flag;

// 			flag = 0;
// 			for (int i = 0; i < subG->AdjList[seq1].size(); i++)
// 			{
// 				if (subG->AdjList[seq1][i] == seq2)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 			{
// 				subG->AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded
// 				subG->eNum++;
// 			}

// 			flag = 0;
// 			for (int i = 0; i < subG->AdjList[seq2].size(); i++)
// 			{
// 				if (subG->AdjList[seq2][i] == seq1)
// 				{
// 					flag = 1;
// 					break;
// 				}
// 			}
// 			if (flag == 0)
// 				subG->AdjList[seq2].push_back(seq1);
// 		}

// 		// M++;
// 		// cout<<M<<endl;
// 	}
// 	unordered_map<int, vector<int>> components;
// 	for (int i : vertices)
// 	{
// 		int root = find(idUFMap[i])->value;
// 		subG->CC[root].insert(i);
// 	}
// }
// bool DataGraph::isCKTruss(vector<int> &vertices, vector<shared_ptr<SkyGroupCand>> &results, int rel, int k, double rho)
// {
// 	DataGraph *relSubG = this->getSubG(vertices);
// 	// TrussDecomposition *kt = new TrussDecomposition(relSubG);
// 	// minimum degree vertex
// 	if (relSubG->eNum < k * (k - 1) / 2)
// 	{
// 		return false;
// 	}
// 	map<Edge, int> support;
// 	int minSup = k - 1;
// 	map<int, shared_ptr<UNode>> idUFMap;
// 	for (int v = 0; v < relSubG->AdjList.size(); v++)
// 	{
// 		if (idUFMap.find(v) == idUFMap.end())
// 		{
// 			shared_ptr<UNode> unode = make_shared<UNode>(v);
// 			makeSet(unode);
// 			idUFMap[v] = unode;
// 		}

// 		vector<int> vNei = relSubG->AdjList[v];
// 		sort(vNei.begin(), vNei.end());
// 		for (int j = 0; j < relSubG->AdjList[v].size(); j++)
// 		{
// 			int u = relSubG->AdjList[v][j];

// 			if (v > u)
// 				continue;
// 			if (idUFMap.find(u) == idUFMap.end())
// 			{
// 				shared_ptr<UNode> unode = make_shared<UNode>(u);
// 				makeSet(unode);
// 				idUFMap[u] = unode;
// 			}
// 			unite(idUFMap[v], idUFMap[u]);
// 			vector<int> uNei = relSubG->AdjList[u];
// 			sort(uNei.begin(), uNei.end());
// 			vector<int> com;
// 			set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
// 			int sup = com.size();
// 			Edge e = make_pair(v, u);
// 			support.emplace(e, sup);
// 			if (sup < minSup)
// 				minSup = sup;
// 		}
// 	}
// 	for (int i = 0; i < relSubG->AdjList.size(); ++i)
// 	{
// 		int root = find(idUFMap[i])->value;
// 		relSubG->CC[root].push_back(i);
// 	}

// 	if (relSubG->CC.size() > 1)
// 	{
// 		if (minSup == k - 2)
// 		{
// 			// find shortest path to satisfy the connected (from edged vertices in a graph to edged vertices in another), then find least vertices to make the added vertices satisfy the k-truss.
// 		}
// 		else if (minSup < k - 2)
// 		{
// 			// find shortest path to satisfy the connected, the find least vertices to make this satisfy the k-truss.
// 		}
// 	}
// 	if (relSubG->CC.size() == 1)
// 	{
// 		if (minSup == k - 2)
// 		{
// 			shared_ptr<SkyGroupCand> p = make_shared<SkyGroupCand>();
// 			p->rel = rel;
// 			p->k = k;
// 			p->rho = rho;
// 			p->vertices = vertices;
// 			results.push_back(p);
// 		}
// 		else if (minSup < k - 2)
// 		{
// 			// select
// 		}
// 	}

// 	// bool isC = isConnected(relSubG->AdjList);
// 	bool isC = false;
// 	if (relSubG->CC.size() == 1)
// 		isC = true;
// 	delete relSubG;
// 	if (minSup == k - 2 && isC)
// 		return true;
// 	else
// 		return false;
// 	// four case:
// 	//  1. ,
// 	//  2. k-truss and disconnected, find shortest path to satisfy the connected (from edged vertices in a graph to edged vertices in another), then find least vertices to make the added vertices satisfy the k-truss.
// 	//  3. not k-truss and connected, find least vertices to make this satisfy the k-truss.
// 	// containing two case: 3.1. k>specified k, which means it is dominated by k, dont need to calculate this k.
// 	//                      3.2. k<specified k, which means it does not dominate k, need to calculate this k.
// 	//  4. not k-truss and disconnected, find shortest path to satisfy the connected, the find least vertices to make this satisfy the k-truss.
// 	// after removing edges whose support is lesser than k-2, the minimum trussness is
// }
// void DataGraph::DFS(int u, vector<int>& disc, vector<int>& low, vector<bool>& visited, int parent, vector<bool>& isCut) {
//     static int time = 0;
//     visited[u] = true;
//     disc[u] = low[u] = ++time;
//     int children = 0;

//     for (int v : AdjList[u]) {
//         if (!visited[v]) {
//             children++;
//             DFS(v, disc, low, visited, u, isCut);
//             low[u] = min(low[u], low[v]);

//             if (parent != -1 && low[v] >= disc[u])
//                 isCut[u] = true;
//         } else if (v != parent) {
//             low[u] = min(low[u], disc[v]);
//         }
//     }

//     if (parent == -1 && children > 1)
//         isCut[u] = true;
// }

// vector<int>& DataGraph::findCutPoints() {
//     vector<int> disc(vNum, -1), low(vNum, -1);
//     vector<bool> visited(vNum, false);
//     vector<bool> isCut(vNum, false);

//     for (int i = 0; i < vNum; i++) {
//         if (!visited[i]) {
//             DFS(i, disc, low, visited, -1, isCut);
//         }
//     }

//     vector<int> cv;
//     for (int i = 0; i < vNum; i++) {
//         if (isCut[i]) {
//             cv.push_back(i);
//         }
//     }
//     return cv;
// }

// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
// discovery time) that can be reached from subtree
// rooted with current vertex
// *st -- >> To store visited edges
// void DataGraph::BCCUtil(set<int> &cv, int ccId, int u, unordered_map<int, int> &disc, unordered_map<int, int> &low, vector<int> &st,
// 						unordered_map<int, int> &parent, vector<vector<int>> &bccs, set<int> &delV)
// {
// 	// A static variable is used for simplicity, we can avoid use
// 	// of static variable by passing a pointer.
// 	static int time = 0;

// 	// Initialize discovery time and low value
// 	disc[u] = low[u] = ++time;
// 	int children = 0;
// 	// store the edge in stack
// 	st.push_back(u);
// 	bool isCut = false;
// 	// Go through all vertices adjacent to this
// 	for (int v : AdjList[u])
// 	{
// 		if (delV.find(v) != delV.end())
// 			continue;
// 		// If v is not visited yet, then recur for it
// 		if (disc[v] == -1)
// 		{
// 			children++;
// 			parent[v] = u;

// 			BCCUtil(cv, ccId, v, disc, low, st, parent, bccs, delV);

// 			// Check if the subtree rooted with 'v' has a
// 			// connection to one of the ancestors of 'u'
// 			// Case 1 -- per Strongly Connected Components Article
// 			low[u] = min(low[u], low[v]);

// 			// If u is an articulation point,
// 			// pop all edges from stack till u -- v
// 			if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u]))
// 			{
// 				isCut = true;
// 				cv.insert(u);
// 				vector<int> bcc;
// 				while (st.back() != u)
// 				{
// 					bcc.push_back(st.back());
// 					st.pop_back();
// 				}
// 				// st.pop_back();
// 				if (!bcc.empty())
// 					bccs.push_back(bcc);
// 			}
// 		}

// 		// Update low value of 'u' only of 'v' is still in stack
// 		// (i.e. it's a back edge, not cross edge).
// 		// Case 2 -- per Strongly Connected Components Article
// 		else if (v != parent[u] && disc[v] < disc[u])
// 		{
// 			low[u] = min(low[u], disc[v]);
// 			// st.push_back(v);
// 			// if (disc[v] < disc[u])
// 			// {
// 			// 	st.push_back(Edge(u, v));
// 			// }
// 		}
// 	}
// 	if (isCut == true)
// 		st.pop_back();
// }
void DataGraph::BCCUtil(set<int> &cv, int u, unordered_map<int, int> &disc, unordered_map<int, int> &low, vector<int> &st,
						unordered_map<int, int> &parent, unordered_map<int, vector<int>> &bccs, set<int> &delV, int &time, int &index)
{
	disc[u] = low[u] = ++time;
	int children = 0;
	st.push_back(u);
	for (int v : AdjList[u])
	{
		if (!disc.count(v) || delV.count(v))
			continue;
		if (disc[v] == -1)
		{
			children++;
			parent[v] = u;

			BCCUtil(cv, v, disc, low, st, parent, bccs, delV, time, index);

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

// The function to do DFS traversal. It uses BCCUtil()
void DataGraph::BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, set<int> &delV, int dv)
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
		if (delV.find(i) != delV.end())
			continue;
		disc[i] = NIL;
		low[i] = NIL;
		parent[i] = NIL;
		// trueSub.push_back(i);
	}

	// if (dv != -1)
	// {
	// 	unordered_set<int> neiDV;
	// 	bool isCore = true;

	// 	for (int u : AdjList[dv])
	// 	{
	// 		if (delV.count(u))
	// 			continue;
	// 		neiDV.insert(u);
	// 	}
	// 	int trigle = 0;
	// 	for (int u : neiDV)
	// 	{
	// 		int neiNum = 0;
	// 		for (int v : AdjList[u])
	// 		{
	// 			if (u > v)
	// 				continue;
	// 			if (!neiDV.count(v))
	// 				continue;
	// 			neiNum++;
	// 			if (neiNum < 2)
	// 			{
	// 				isCore = false;
	// 				break;
	// 			}
	// 			for (int w : AdjList[v])
	// 			{
	// 				if (neiDV.count(w) && count(AdjList[u].begin(), AdjList[u].end(), w))
	// 					trigle++;
	// 			}
	// 		}
	// 		if (!isCore)
	// 			break;
	// 	}

	// 	// if (isCore == true && trigle >= 3 * (neiDV.size() - 2))
	// 	// {
	// 	// 	bccs[index] = trueSub;
	// 	// 	cout << "fifth pruning." << endl;
	// 	// 	return;
	// 	// }
	// }

	for (int i : subG)
	{
		auto t = disc.find(i);
		if (t != disc.end() && t->second == NIL)
			BCCUtil(cv, i, disc, low, st, parent, bccs, delV, time, index);

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
void DataGraph::support(vector<Edge> &kSupEs, int &minSup, map<Edge, int> &sups, int ccId, set<int> &delV, int k)
{
	for (int v : CC[ccId])
	{
		if (delV.count(v))
			continue;

		set<int> vNei;
		for (int u : AdjList[v])
		{
			if (!delV.count(u))
				vNei.insert(u);
		}
		for (int u : vNei)
		{
			if (v > u)
				continue;
			set<int> uNei;
			for (int i : AdjList[u])
			{
				if (!delV.count(i))
					uNei.insert(i);
			}
			set<int> com;
			set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), inserter(com, com.begin()));
			int sup = com.size();
			Edge e = make_pair(v, u);
			sups.emplace(e, sup);
			minSup = min(sup, minSup);
			if (sup == k - 2)
				kSupEs.push_back(e);
			// if(!sup2edge.count(sup)){
			// 	sup2edge.emplace(sup, vector<Edge>());
			// }
			// sup2edge[sup].push_back(e);
		}
	}
}
void DataGraph::support(int k, int ccId, unordered_map<int, unordered_map<int, int>> &e2sup, queue<Edge> &supLesserKE, unordered_set<int> &delV)
{
	for (int v : CC[ccId])
	{

		if (delV.count(v))
			continue;

		set<int> vNei;
		for (int u : AdjList[v])
		{
			if (!delV.count(u))
				vNei.insert(u);
		}
		for (int j = 0; j < AdjList[v].size(); j++)
		{
			int u = AdjList[v][j];
			if (v > u)
				continue;

			set<int> uNei;
			for (int i : AdjList[u])
			{
				if (!delV.count(i))
					uNei.insert(i);
			}
			vector<int> com;
			set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
			int s = com.size();

			auto it1 = e2sup.find(v);
			if (it1 != e2sup.end())
			{
				it1->second.emplace(u, s);
			}
			else
				e2sup.emplace(v, unordered_map<int, int>{{u, s}});
			auto it2 = e2sup.find(u);
			if (it2 != e2sup.end())
			{
				it2->second.emplace(v, s);
			}
			else
				e2sup.emplace(u, unordered_map<int, int>{{v, s}});
			if (s < k - 2)
			{
				Edge e = make_pair(v, u);
				supLesserKE.push(e);
			}
		}
	}
}
void DataGraph::support(vector<Edge> &kSupEs, int &minSup, unordered_map<int, unordered_map<int, int>> &e2sup, int ccId, int k)
{

	for (int v : CC[ccId])
	{
		sort(AdjList[v].begin(), AdjList[v].end());
		for (int u : AdjList[v])
		{
			if (v > u)
				continue;
			set<int> uNei;
			sort(AdjList[u].begin(), AdjList[u].end());
			set<int> com;
			set_intersection(AdjList[v].begin(), AdjList[v].end(), AdjList[u].begin(), AdjList[u].end(), inserter(com, com.begin()));
			int sup = com.size();
			auto it1 = e2sup.find(v);
			if (it1 != e2sup.end())
			{
				it1->second.emplace(u, sup);
			}
			else
				e2sup.emplace(v, unordered_map<int, int>{{u, sup}});
			auto it2 = e2sup.find(u);
			if (it2 != e2sup.end())
			{
				it2->second.emplace(v, sup);
			}
			else
				e2sup.emplace(u, unordered_map<int, int>{{v, sup}});
			minSup = min(sup, minSup);
			if (sup == k - 2)
				kSupEs.push_back(make_pair(v, u));
		}
	}
}
void DataGraph::support(map<Edge, int> &sups, map<int, vector<Edge>> &sup2edge, set<int> &subG) // subG can be unordered_set
{
	for (int v : subG)
	{
		set<int> vNei;
		for (int u : AdjList[v])
		{
			if (subG.count(u))
				vNei.insert(u);
		}
		for (int u : vNei)
		{
			if (v > u)
				continue;
			set<int> uNei;
			for (int i : AdjList[u])
			{
				if (subG.count(u))
					uNei.insert(i);
			}
			set<int> com;
			set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), inserter(com, com.begin()));
			int sup = com.size();
			Edge e = make_pair(v, u);
			sups.emplace(e, sup);
			if (!sup2edge.count(sup))
			{
				sup2edge.emplace(sup, vector<Edge>{e});
			}
			else
				sup2edge[sup].push_back(e);
		}
	}
}
void DataGraph::support(unordered_map<int, unordered_map<int, unordered_map<int, int>>> &e2sups, unordered_map<int, vector<Edge>> &sup2edge)
{
	for (auto cc : CC)
	{
		int ccId = cc.first;
		unordered_map<int, unordered_map<int, int>> e2sup;
		for (int v : cc.second)
		{
			vector<int> vNei = AdjList[v];
			sort(vNei.begin(), vNei.end());
			for (int j = 0; j < AdjList[v].size(); j++)
			{
				int u = AdjList[v][j];
				if (v > u)
					continue;

				vector<int> uNei = AdjList[u];
				sort(uNei.begin(), uNei.end());
				vector<int> com;
				set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
				int s = com.size();
				auto it1 = e2sup.find(v);
				if (it1 != e2sup.end())
				{
					it1->second.emplace(u, s);
				}
				else
					e2sup.emplace(v, unordered_map<int, int>{{u, s}});
				auto it2 = e2sup.find(u);
				if (it2 != e2sup.end())
				{
					it2->second.emplace(v, s);
				}
				else
					e2sup.emplace(u, unordered_map<int, int>{{v, s}});
				Edge e = make_pair(v, u);
				auto it3 = sup2edge.find(s);
				if (it3 == sup2edge.end())
				{
					sup2edge.emplace(s, vector<Edge>{e});
				}
				else
					it3->second.push_back(e);
			}
		}
		e2sups.emplace(ccId, e2sup);
	}
	// for (int v = 0; v < AdjList.size(); v++)
	// {
	//     vector<int> vNei = AdjList[v];
	//     sort(vNei.begin(), vNei.end());
	//     for (int j = 0; j < AdjList[v].size(); j++)
	//     {
	//         int u = AdjList[v][j];
	//         if (v > u)
	//             continue;
	//         vector<int> uNei = AdjList[u];
	//         sort(uNei.begin(), uNei.end());
	//         vector<int> com;
	//         set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
	//         int s = com.size();
	//         auto it1 = e2sups.find(v);
	//         if (it1 != e2sups.end())
	//         {
	//             it1->second.emplace(u, s);
	//         }
	//         else
	//             e2sups.emplace(v, unordered_map<int, int>{{u, s}});
	// 		auto it2 = e2sups.find(u);
	//         if (it2 != e2sups.end())
	//         {
	//             it2->second.emplace(v, s);
	//         }
	//         else
	//             e2sups.emplace(u, unordered_map<int, int>{{v, s}});
	//     }
	// }
}
void DataGraph::support(unordered_map<int, unordered_map<int, int>> &e2sup, set<int> &sub)
{
	unordered_map<int, set<int>> v2Nei;
	for (int v : sub)
	{
		vector<int> vWholeNei = AdjList[v];
		sort(vWholeNei.begin(), vWholeNei.end());
		set<int> vNei;
		set_intersection(vWholeNei.begin(), vWholeNei.end(), sub.begin(), sub.end(), inserter(vNei, vNei.begin()));
		v2Nei.emplace(v, vNei);
	}
	for (int v : sub)
	{
		set<int> &vNei = v2Nei[v];
		for (int u : vNei)
		{
			if (v > u)
				continue;

			set<int> &uNei = v2Nei[u];
			vector<int> com;
			set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
			int s = com.size();
			e2sup[v][u] = s;
			e2sup[u][v] = s;
		}
	}
}
