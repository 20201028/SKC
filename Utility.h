#pragma once
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
#include "DataGraph.h"
void dfs(int node, unordered_set<int> &visited, unordered_map<int, unordered_map<int, int>> &G);
void firstMaxRhoLesser(unordered_map<int, unordered_map<int, int>> &shrinkG, int k, int ccId, DataGraph *G, unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho);
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, set<int> &delV, unordered_map<int, unordered_map<int, int>> &edge2sup);
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, set<int> &delV, int dv, DataGraph *relSubG);
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, set<int> &subG, unordered_map<int, unordered_map<int, int>> &g);
void BCC(set<int> &cv, unordered_map<int, vector<int>> &bccs, vector<int> &subG, unordered_map<int, unordered_map<int, int>> &g);
bool checkInd(int v, int k, unordered_map<int, unordered_map<int, int>> &edge2sup);
bool checkInd(int v, DataGraph *g, set<int> &removeV, int k, unordered_map<int, unordered_map<int, int>> &edge2sup);
bool checkInd(int v, DataGraph *g, set<int> &removeV, int k, map<Edge, int> &edge2sup, vector<Edge> &kSupEs);
bool checkInd(int v, int k, unordered_map<int, unordered_map<int, int>> &edge2sup, vector<Edge> &kSupEs);
bool dominates(const SkyGroupCand &a, const SkyGroupCand &b);
void skyline(vector<SkyGroupCand> &result, const vector<SkyGroupCand> &points);
bool firstMaxRhoLesser(unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho);
bool secondMaxRhoLesser(unordered_map<int, double> &vert2weightCopy, unordered_map<int, unordered_map<int, double>> &simG, double maxRho, int delV);
void updateBcc(set<int> &cv, unordered_map<int, vector<int>> &bccId2verts, set<int> &cv1, unordered_map<int, vector<int>> &bccId2verts1, unordered_map<int, vector<int>> &vert2bccId, int bccId, int &bccIndex);