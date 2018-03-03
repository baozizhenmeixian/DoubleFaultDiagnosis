#pragma once
/* 2013.8 By ����ά
* �������ڴ��һЩ���õķ���
*/
#include <string>
#include <vector>
#include <sstream>
using namespace std;

class Utility
{
public:
	Utility(void);
	~Utility(void);
	string getLiteralNegation(string literal);
	string sortTerm(string term);
	vector<string> uniqueTermVector(vector<string> termVector, string term);
	bool evaluate(string exp,int testData[],string model);
	bool evaluate(string exp,int testNum);
	bool evaluate(string exp, int testNum, string model);
	bool evaluate(string exp,string testString);
	bool evaluateMutant(string oriExp,string mutant,string testString);
	bool evaluateDNF(string exp,int testData[],string model);
	bool evaluateDNF(string exp,int testNum, string model);
	bool evaluateDNF(string exp,string testString,string model);
	vector<string> extractCNFstyle(string exp);
	vector<string> extractDNFstyle(string exp);
	string generateModel(string exp);
	void intToIntArray(int intNum,int intArray[],int arrayLength);
	string intTostring(int intNum,int length);
	bool isCNFtermFalse(string term,int testData,string model);
	bool isCNFtermFalse(string term,int testData[],string model);
	bool isDNFtermTrue(string term,int testData,string model);
	bool isDNFtermTrue(string term,int testData[],string model);
	string sortVariables(string variables);
	string sortDNF(string dnf);
	string sortCNF(string cnf);
	bool isTwoPointsAdjacent(int p1[],int p2[],int len);
	bool isTwoPointsAdjacent(int p1,int p2,int len);
	vector<vector<int>> getConnectedRegion(string exp, bool trueOrFalsePoint);
	vector<int> getTermsOpCount(string exp,bool DNF1CNF0);
	bool isExpressionsEquivalent(string exp1, string exp2);

	//�Ա���ʽ����Ԥ����������ŷֵ�terms�У�vars��¼��ĸ�Ƿ����
	void mutantsPreproccess(string exp, vector<string> &terms, vector<vector<string>> &literals_in_terms, bool (&vars)[26]);
	//�滻����ʽ�е�����,expΪԭ����ʽ��a,bΪԭ����ʽ����Ҫ���滻��term���±�
	string replaceDoubleTerms(string exp, int a,int b,string term1,string term2);


	string replaceDoubleTerms(string exp, vector<string> terms, int a, int b, string term1, string term2);

	int getTermPosition(vector<string> terms, string term);
};
