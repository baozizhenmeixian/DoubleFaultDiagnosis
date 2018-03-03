#pragma once
#include "CommonFuncs.h"
using namespace std;
class CDLRFHandler
{
public:
	CDLRFHandler(void);
	~CDLRFHandler(void);

	CCommonFuncs m_commonFuncs;

	//����LRF�㷨
	string run(string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, bool mutantsNeed,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests, vector<int> allOTPs, vector<int> allUTPs);
	string run(
		string dimensionValuesStr,
	bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<int> allOTPs, 
	vector<int> allUTPs
	);
	string run2(
		string dimensionValuesStr,
		bool vOdd,
		vector<vector<CCell*>>& matrixCells,
		vector<int> allOTPs,
		vector<int> allUTPs
		);
	string runWithParameter(int param, vector<int>&literalNum, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allUTPs);
	string runWithParameter2(int param, vector<int>&literalNum, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allUTPs);
	string runWithParameter(int param, vector<int>&literalNum, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, bool mutantsNeed,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests, vector<int> allOTPs, vector<int> allUTPs);
	//��ʾϵ������
	string getShowStr(vector<vector<CCell*>> matrixCells);
};

