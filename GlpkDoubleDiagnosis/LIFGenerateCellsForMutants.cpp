#pragma once
#include "stdafx.h"
#include "LIFGenerateCellsForMutants.h"
using namespace std;
string LIFGenerateCellsForMutants(
	CCommonFuncs commonFuncs,
	string dimensionValuesStr,
	bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests,
	vector<int> allOTPs)
{
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;
	//��ȡ��άֵ
	commonFuncs.getDimensionValue(dimensionValuesStr, vOdd, dimensionValues, biggestDimensionValueAllowed);
	int dimensionNum = dimensionValues.size();
	vector<vector<CCell*>> shrinkCells;
	//��ÿһά������С����ȡ��С���ֵĵ�
	for (int i = 0; i<dimensionNum; i++)
	{
		vector<CCell*> shrinkCell;
		vector<vector<int>> tempShrinkDimensionValue;
		vector<int> shrinkDimensionValue;
		int ithDimensionNum = dimensionValues.at(i).size();
		switch (ithDimensionNum)
		{
		case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
			for (int m = 1; m<5; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(commonFuncs.boundValue(m, 4));
				shrinkDimensionValue.push_back(commonFuncs.boundValue(m + 1, 4));
				//tempShrinkDimensionValue  : �̺���ǰ�����ȱ�ݸ񼯺�
				commonFuncs.getTempDimensionValue(dimensionValues, i, shrinkDimensionValue, tempShrinkDimensionValue);
				commonFuncs.ListAllCells(0, tempShrinkDimensionValue, shrinkCell);
				commonFuncs.kickOutCells(shrinkCell, allOTPs);
				if (shrinkCell.size()>0){
					shrinkCells.push_back(shrinkCell);
				}
				// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(commonFuncs.boundValue(m + 2, 4));
				shrinkDimensionValue.push_back(commonFuncs.boundValue(m + 3, 4));
				vector<vector<int>> tempMutantDimensionValue;
				commonFuncs.getTempDimensionValue(dimensionValues, i, shrinkDimensionValue, tempMutantDimensionValue);
				mutantsDimensionValue.push_back(tempMutantDimensionValue);

				//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
				corresponding_tests.push_back(shrinkCell);
			}
			break;
		case 2: //�����������Ϊ1
			for (int m = 0; m<2; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(commonFuncs.boundValue(dimensionValues.at(i).at(m), 4));
				commonFuncs.getTempDimensionValue(dimensionValues, i, shrinkDimensionValue, tempShrinkDimensionValue);
				commonFuncs.ListAllCells(0, tempShrinkDimensionValue, shrinkCell);
				commonFuncs.kickOutCells(shrinkCell, allOTPs);//Shrink�ĵ㲻����OTP
				if (shrinkCell.size()>0){
					shrinkCells.push_back(shrinkCell);
				}

				// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(commonFuncs.boundValue(dimensionValues.at(i).at((m + 1) % 2), 4));
				vector<vector<int>> tempMutantDimensionValue;
				commonFuncs.getTempDimensionValue(dimensionValues, i, shrinkDimensionValue, tempMutantDimensionValue);
				mutantsDimensionValue.push_back(tempMutantDimensionValue);

				//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
				corresponding_tests.push_back(shrinkCell);
			}
			break;
		case 1://������������ٽ�������
			break;
		}
	}

	//�Ѿ���ȡ���е�shrink cells�������п��ܵ�shrink���һ������һ�б�ʾһ��shrink
	for (int i = 0; i<shrinkCells.size(); i++)
	{
		vector<CCell*> matrixRow;
		for (int j = 0; j<shrinkCells.at(i).size(); j++)
			matrixRow.push_back(shrinkCells.at(i).at(j));
		matrixCells.push_back(matrixRow);
	}

	return "";
}

