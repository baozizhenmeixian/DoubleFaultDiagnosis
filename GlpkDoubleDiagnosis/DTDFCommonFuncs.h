#pragma once
#include "CommonFuncs.h"
using namespace std;
class DTDFCommonFuncs
{
public:
	DTDFCommonFuncs(void);
	~DTDFCommonFuncs(void);

	CCommonFuncs m_commonFuncs;

	// ����LIFxLIF�㷨����Լ������
	string generateMatrixofOTP(string dimensionValuesStr1, string dimensionValuesStr2, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> all3OTPs, bool &hasequal, bool vOdd);
	string generateMatrixofRFP(string dimensionValuesStr1, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allFPs, bool &hasequal, bool vOdd);
};