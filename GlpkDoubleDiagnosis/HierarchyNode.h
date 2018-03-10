#pragma once
#include <string>
using namespace std;

class HierarchyNode
{
public:
	HierarchyNode(string);
	~HierarchyNode(void);

	void setValue(string v){ value = v; };
	string getValue()const{ return value; };

	void setNextCondition(vector<int> v){ nextConditon = v; };
	vector<int> getNextCOndition(){ return nextConditon; };

	void setIsOut(bool v){ isOut = v; };
	int getIsOut()const{ return isOut; };

	void setNext(vector<HierarchyNode> v){ next = v; };
	vector<HierarchyNode> getNext(){ return next; };





private:
	string value;
	vector<HierarchyNode> next;
	vector<int> nextConditon;
	bool isOut;
};