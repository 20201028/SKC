#pragma once
#include <memory>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#define NIL -1
#define denStyle 0.0
// #include <limits.h>
using namespace std;
using Item = int;
// using Transaction = set<Item>;
// using TransformedPrefixPath = std::pair<std::vector<Item>, uint64_t>;
// using TransformedPrefixPath = pair<Transaction, uint64_t>;
// using Pattern = std::pair<std::set<Item>, uint64_t>;
// using TransSet = vector<int>;
// using TransformedPrefixPath = std::pair<std::vector<Item>, TransSet>;
// using Pattern = pair<set<Item>, TransSet>;
using Edge = pair<int, int>;
struct DEdge {
    int to;     // 目标节点
    double weight; // 边的权重
};
using EdgesSet = vector<Edge>;
// using Clique = vector<int>;
using Frac = pair<int, int>;
// #ifndef MIN
//   #define MIN(a,b) ((a)<(b)?(a):(b))
// #endif

// #ifndef MAX
//   #define MAX(a,b) ((a)>(b)?(a):(b))
// #endif

extern "C"
{
    double densest_subgraph(int64_t n, int64_t m, vector<int> &ei, vector<int> &ej, vector<double> &w, vector<int> &output);
};
// int gcd(int num, int deno);
// struct Weight
// {
//     int numerator;
//     int denominator;
//     Weight()
//     {

//     }
//     Weight(int num, int deno)
//     {
//         int c = gcd(num, deno);
//         numerator = num / c;
//         denominator = deno / c;
//     }
//     bool operator < (const Weight& b) const
//     {
//         return numerator * b.denominator - denominator * b.numerator < 0;
//     }
//     bool operator > (const Weight& b) const
//     {
//         return numerator * b.denominator - denominator * b.numerator > 0;
//     }
//     bool operator == (const Weight& b) const
//     {
//         return numerator == b.numerator && denominator == b.denominator;
//     }
//     Weight operator * (const Weight& b) const
//     {
//         return Weight(numerator * b.numerator, denominator * b.denominator);
//     }
//     Weight operator * (const int& b) const
//     {
//         return Weight(numerator * b, denominator);
//     }
//     Weight Dist() const // ��1 - a/b
//     {
//         return Weight(denominator - numerator, denominator);
//     }
// };
struct SkyGroupCand{
    int rel;
    int k;
    double rho;
    set<int> vertices;
};
struct UNode
{
	int value;
	shared_ptr<UNode> parent = nullptr;
	int rank = -1;

	// Constructor
	UNode(int value) : value(value) {}
};

void makeSet(shared_ptr<UNode> &x);

shared_ptr<UNode> find(shared_ptr<UNode> &x);

void unite1(shared_ptr<UNode> &x, shared_ptr<UNode> &y);



