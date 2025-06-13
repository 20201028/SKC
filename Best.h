#include "DataGraph.h"
#include "KTruss.h"
#include "Utility.h"
// 定义一个结构表示图中的边


// 定义一个优先队列中的元素
struct Node {
    int vertex;   // 当前节点
    double distance; // 到该节点的距离

    // 定义优先队列的比较规则（距离小的优先）
    bool operator>(const Node& other) const {
        return distance > other.distance;
    }
};
// double dijkstra(int source, int target, unordered_map<int, unordered_map<int, double>>& graph, vector<int>& predecessor);
// void printPath(int source, int target, const vector<int>& predecessor);
void Best(map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
void Best(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);
void Best1(int k, map<int, vector<Edge>, greater<int>> &rel2edge, vector<SkyGroupCand> &result, map<int, unordered_map<int, double>> &simE, unordered_map<int, set<int>> &rel2verts);