#pragma once
#include <vector>
#include "biddy.h"
using namespace std;

class BDD
{
public:
	//generate instance
	BDD(string expression);
	BDD(BDD& bdd);
	~BDD(void);

	//��ֵ�ӿ�
	int evaluate(string testStr);
	int evaluate(int testInt);
	int evaluate(std::vector<int> testIntArr);

	//��Ϊ������ֵ
	int evaluateAsMutant(string testStr, string oriExpVars);
	int evaluateAsMutant(int testInt, string oriExpVars);
	int evaluateAsMutant(std::vector<int> testIntArr, string oriExpVars);
	
	//�ж������������ʽ�Ƿ�ȼ�
	bool isEquivalent(BDD& bdd);

	//
	string getPrintString();

	//
	string getExpression()const{return m_expression;};
	string getVariables()const{return m_variables;};
	Biddy_Edge getBiddyEdge()const{return m_BDD;};

private:
	static int _instanceCount;
	int generateBDD();

	//member variables
	string m_expression;
	string m_variables;
	std::vector<Biddy_Edge> m_variableNodes;
	Biddy_Edge m_BDD;
};




//#define DEBUG_BDD

//*************//unit test
#ifdef DEBUG_BDD

class BDD_UnitTest
{
public:
	BDD_UnitTest();
	~BDD_UnitTest();

	void run();

	void testEvaluate();
	void testEvaluateMutant();
	void testIsEquivalent();
private:

};

#endif
