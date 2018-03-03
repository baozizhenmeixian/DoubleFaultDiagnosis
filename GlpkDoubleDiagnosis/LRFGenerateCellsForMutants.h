#pragma once
#include "CommonFuncs.h"
using namespace std;
string LRFGenerateCellsForMutants(
	CCommonFuncs commonFuncs,
	string dimensionValuesStr,
	bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests,
	vector<int> allOTPs,
	vector<int> allUTPs
);
string runWithParameter(
	CCommonFuncs commonFuncs,
	int param,
	vector<int>&literalNum, 
	string dimensionValuesStr, 
	bool vOdd, 
	vector<vector<CCell*>>& matrixCells,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests,
	vector<int> allOTPs, vector<int> allUTPs
);

