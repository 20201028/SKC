#pragma once
#include "Define.h"
#include <map>
#include <unordered_map>
#include "KTruss.h"
#include "Utility.h"

void Improve(map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, 
map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
void Improve(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, 
map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
void Improve1(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, 
map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
