#pragma once
#include "Define.h"
#include <map>
#include <unordered_map>

#include "Utility.h"
#include "KTruss.h"
void Basic(map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
void Basic(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
void Basic1(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);