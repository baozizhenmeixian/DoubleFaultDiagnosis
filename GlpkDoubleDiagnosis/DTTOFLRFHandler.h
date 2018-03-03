#pragma once
#include "CommonFuncs.h"
using namespace std;
class CDTTOFLRFHandler
{
public:
	CDTTOFLRFHandler(void);
	~CDTTOFLRFHandler(void);
	CCommonFuncs m_commonFuncs;

	//����LIFxLRF�㷨
	string run(string dimensionValuesStr1, string dimensionValuesStr2,
		bool vOdd, 
		vector<vector<CCell*>>& matrixCells, 
		vector<int> allOTPs,vector<int> allUTPs);

	string run2(string dimensionValuesStr1, string dimensionValuesStr2,
		bool vOdd,
		vector<vector<CCell*>>& matrixCells,
		vector<int> allOTPs, vector<int> allUTPs);
};

