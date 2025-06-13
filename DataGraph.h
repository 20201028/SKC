#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <map>
#include <iostream>
#include <set>
#include <list>
#include <queue>
#include "Define.h"

using namespace std;

// struct SimG
// {
// 	int n, m;
// 	vector<int> ei, ej;
// 	vector<double> w;
// };
class DataGraph
{
public:
	/*
	 * dimension 1: vertex ID := iterator sequence
	 * dimension 2: the list of adjacent verticess
	 */
	int vNum,eNum;
	vector<vector<int>> AdjList;
	// vector<int> verSort;
	unordered_map<int, int> id2seq;		 // map the vertex from ordinary dataset id to AdjList sequence. id of a node to line index of this node in adjacency list
	unordered_map<int, int> seq2id;		 // map the sequence back to the id for presentation convenience
	unordered_map<int, set<int>> id2att; // seq and original attribute
	// vector<unordered_map<int, double>> SimAdjList;
	unordered_map<int, int> seq2rel;
	map<int, set<int>> rel2id;
	unordered_map<int, set<int>> CC;
	unordered_map<int,map<int, unordered_map<int, double>>> SimG;
	map<int, shared_ptr<UNode>> idUFMap;
	// unordered_map<int, int> id2att;//index:seq,element:attribute

	DataGraph(string &attGraph, int &N, int &M);
	
	// DataGraph(string &attGraph, int &N, int &M);
	DataGraph(string &realGraph, string &imagGraph, int &N, int &M);
	DataGraph(unordered_map<int, int> &, DataGraph *);
	DataGraph(vector<pair<int, int>> &);
	DataGraph();
	// DataGraph *getSubG(set<pair<int, int>> &);
	DataGraph *getSubG(vector<int> &vertices);
	void getSubG(unordered_map<int,set<int>> &cc, unordered_map<int, unordered_map<int, int>> &edge2sup, set<int> &vertices);
	unordered_map<int,map<int, unordered_map<int, double>>> &getSimGs(unordered_map<int, set<int>> &);
	void getSimGs(double simThreshold, map<int, unordered_map<int, double>> &, unordered_map<int, set<int>> &, int);
	// bool isKTruss(vector<int> &vertices, int k);
	void getCC(unordered_map<int, set<int>> &ccDS, set<int> &maxRV);
	bool isCKTruss(vector<int> &vertices, vector<shared_ptr<SkyGroupCand>> &results, int k, int rel, double rho);
	void addEdges(vector<Edge> &, vector<pair<int, int>> &);
	void addEdgesAndMaintainC(vector<pair<int, int>> &edges);
	void addEdge(Edge &);
	void addEdge(int src, int dst);
	void addEdgeNoMatinC(int src, int dst);
	void addEdgeAndMainConnect(int src, int dst, map<int, unordered_map<int, double>> &simE, unordered_map<int, unordered_map<int, double>> &simG);
	void addEdgeAndMainConnect(vector<Edge> &newEs, int src, int dst, map<int, unordered_map<int, double>> &simE, unordered_map<int, unordered_map<int, double>> &simG);
	void unite(shared_ptr<UNode> &x, shared_ptr<UNode> &y);
	void DFS(int u, vector<int>& disc, vector<int>& low, vector<bool>& visited, int parent, vector<bool>& isCut);
    vector<int>& findCutPoints();
	void BCC(set<int> &cv, unordered_map<int, vector<int>> &, vector<int> &subG, set<int> &removeV, int dv);
	void support(int k, int ccId, unordered_map<int, unordered_map<int, int>> &e2sup, queue<Edge> &supLesserKE, unordered_set<int> &delV);
	void support(vector<Edge> &kSupEs, int &minSup, map<Edge, int> &sups, int ccId, set<int> &removeV, int k);
	void support(vector<Edge> &kSupEs, int &minSup, unordered_map<int, unordered_map<int, int>> &sups, int ccId, int k);
	void support(map<Edge, int> &sups, map<int, vector<Edge>> &sup2edge, set<int> &subG);
	void support(unordered_map<int,unordered_map<int, unordered_map<int, int>>> &e2sups, unordered_map<int,vector<Edge>> &sup2edge);
	void support(unordered_map<int, unordered_map<int, int>> &e2sup, set<int> &sub);
	void BCCUtil(set<int> &cv, int u, unordered_map<int,int> &disc, unordered_map<int,int> &low, vector<int> &st,
					unordered_map<int,int> &parent, unordered_map<int, vector<int>> &bccs, set<int> &delV, int &time, int &index);
	// void BCC(int ccId);
	// void BCCUtil(int ccId, int u, unordered_map<int,int> &disc, unordered_map<int,int> &low, vector<Edge> &st,
	// 				unordered_map<int,int> &parent);
	// void addVert(int &, DataGraph*);
	// bool compVertex(int i, int j);
};
