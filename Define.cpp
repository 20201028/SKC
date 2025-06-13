#include "Define.h"
void makeSet(shared_ptr<UNode> &x)
{
	x->parent = x;
	x->rank = 0;
}
shared_ptr<UNode> find(shared_ptr<UNode> &x)
{
	if (x->parent != x)
	{
		x->parent = find(x->parent); // make the parent of x point to the root
	}
	return x->parent;
}
void unite1(shared_ptr<UNode> &x, shared_ptr<UNode> &y) // union
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
	}
	else if (xRoot->rank > yRoot->rank)
	{
		yRoot->parent = xRoot;
	}
	else
	{
		yRoot->parent = xRoot;
		xRoot->rank = xRoot->rank + 1;
	}
}