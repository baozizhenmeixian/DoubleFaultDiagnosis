#pragma once
#include "StdAfx.h"
#include "LIFGenerateCellsForMutants.h"
//����LRF�㷨  paramΪ1��ֻ����validShrink����  paramΪ2��ֻ����expand���� ����ֵ����������
string runWithParameter(
	CCommonFuncs commonFuncs,
	int param,
	vector<int>&literalNum,
	string dimensionValuesStr,
	bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests,
	vector<int> allOTPs,
	vector<int> allUTPs
	)
{
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;
	int allShrinksCount = 0;
	int validShrinksCount = 0;

	// if mutantsNeed ,how to identify a equivalent LIF mutant 
	// Then if yes, generate a TIF_LRF_mutants 
	commonFuncs.getDimensionValue(dimensionValuesStr, vOdd, dimensionValues, biggestDimensionValueAllowed);

	int dimensionNum = dimensionValues.size();
	//���ÿһά����Ϊԭ����1/2����ѡһά����Ϊԭ��������

	for (int i = 0; i<dimensionNum; i++)
	{
		//�ֱ��ȡʣ�µ�һ���ȥ����һ��
		vector<vector<int>> remainedDimensionValues;
		vector<vector<int>> removedDimensionValues;
		vector<int> remainedTemp;
		vector<int> removedTemp;

		int ithDimensionValueNum = dimensionValues.at(i).size();
		//�о���һά���ܳ��ֵ������������¼����
		switch (ithDimensionValueNum)
		{
		case 4://��һά��1,2,3,4�����
			// Let XY represent the i-th dimension
			/*   to remove      will remain   meaning(insert a ...)
			*    1 , 2          3 , 4          insert  X
			*    2 , 3          1 , 4          insert  !Y
			*    3 , 4          1 , 2		   insert  !X
			*    4 , 1          2 , 3		   insert  Y
			*/
			for (int m = 1; m<5; m++)
			{
				remainedTemp.clear();
				removedTemp.clear();
				remainedTemp.push_back(commonFuncs.boundValue(m, 4));
				remainedTemp.push_back(commonFuncs.boundValue(m + 1, 4));
				removedTemp.push_back(commonFuncs.boundValue(m + 2, 4));
				removedTemp.push_back(commonFuncs.boundValue(m + 3, 4));
				remainedDimensionValues.push_back(remainedTemp);
				removedDimensionValues.push_back(removedTemp);
			}
			break;
		case 2://��һά��1,2��2,3��3,4��4,1�����
			remainedTemp.clear();
			removedTemp.clear();
			remainedTemp.push_back(dimensionValues.at(i).at(0));
			removedTemp.push_back(dimensionValues.at(i).at(1));
			remainedDimensionValues.push_back(remainedTemp);
			removedDimensionValues.push_back(removedTemp);

			remainedTemp.clear();
			removedTemp.clear();
			remainedTemp.push_back(dimensionValues.at(i).at(1));
			removedTemp.push_back(dimensionValues.at(i).at(0));
			remainedDimensionValues.push_back(remainedTemp);
			removedDimensionValues.push_back(removedTemp);

			break;
		case 1://��һάֻ��һ��ֵ�����ܽ�������
			break;
		}

		//���ÿһ�ֿ��ܳ��ֵ�����������µ�����Ϊ��������ȡԼ����ȥ����һ��������ĵ㣩
		for (int j = 0; j<removedDimensionValues.size(); j++)
		{
			//��ȡȥ������һ��ĵ㣬�����һ��ȥ���ж�������
			vector<CCell*> removedCells;
			vector<vector<int>> tempDimensionValue;

			//��ȡ��ȥ����һ��ĵ�
			commonFuncs.getTempDimensionValue(dimensionValues, i, removedDimensionValues.at(j), tempDimensionValue);
			commonFuncs.ListAllCells(0, tempDimensionValue, removedCells);
			commonFuncs.kickOutCells(removedCells, allOTPs);//ȥ����һ���в���ѡ��OTP
			allShrinksCount++;
			bool hasEquivalentLifOccurred = false;
			/*
			*
			*ȥ��OTP֮������жϣ�
			*1)���removedCells�л��е㣬��˵����shrink��һ����Чshrink������������ĵ㡣
			*2)���removedCellsΪ�ռ�������Ҫ����������ĵ�
			*/

			if (removedCells.size() > 0){//��Ч�����������¼
				hasEquivalentLifOccurred = true;
				//����Ӧ���ŵ���Ҫ��¼Լ������
			}
			else{//��Ч��������Ҫ��¼
				hasEquivalentLifOccurred = false;
			}

			if (!hasEquivalentLifOccurred){
				//��ȡʣ�µ�һ���dimensionֵ
				vector<vector<int>> tempRemainedDimensionValue;
				commonFuncs.getTempDimensionValue(dimensionValues, i, remainedDimensionValues.at(j), tempRemainedDimensionValue);

				//�������ȥ�����г����п��ܵ�����
				vector<vector<CCell*>> expandCells;
				int literalN = 0;//
				for (int k = 0; k<dimensionNum; k++)//�����ж�k==i�����
				{
					vector<CCell*> expandCell;
					vector<vector<int>> tempExpandDimensionValue;
					vector<int> expandDimensionValue;
					int biggestValue = biggestDimensionValueAllowed.at(k);//��ά��������ֵ��2��4��
					int kthRemainedDimensionValueNum = tempRemainedDimensionValue.at(k).size();
					switch (kthRemainedDimensionValueNum)
					{
					case 4://�����������һά������������
						break;
					case 2://����������������ʣ�µ�����������ԭ��Ϊ1,2�����䲿�־�Ϊ3,4
						if (k != i&&biggestValue == 4)//1.k==iʱexpand֮���shrink֮ǰ��һ����;2.�����ά�����ȡֵֻ��2��Ҳ���ܽ���expand
						{
							/*Let XY represent the k-th dimension
							*    init          expanded        meaning(omit a ...)
							*    1 , 2          3 , 4          omit  X
							*    2 , 3          1 , 4          omit  !Y
							*    3 , 4          1 , 2		   omit  !X
							*    4 , 1          2 , 3		   omit  Y
							*/

							int num1 = tempRemainedDimensionValue.at(k).at(0);
							int num2 = tempRemainedDimensionValue.at(k).at(1);
							if ((num1 == 1 && num2 == 4) || (num1 == 4 && num2 == 1))
							{
								expandDimensionValue.push_back(2);
								expandDimensionValue.push_back(3);
							}
							else
							{
								int minNum = num1<num2 ? num1 : num2;
								expandDimensionValue.push_back(commonFuncs.boundValue(minNum + 2, biggestValue));
								expandDimensionValue.push_back(commonFuncs.boundValue(minNum + 3, biggestValue));
							}
							commonFuncs.getTempDimensionValue(tempRemainedDimensionValue, k, expandDimensionValue, tempExpandDimensionValue);

							commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
							commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
							if (removedCells.size() + expandCell.size()>0){
								expandCells.push_back(expandCell);
								literalNum.push_back(literalN);
							}
							literalN++;

							vector<CCell*> matrixRow;
							// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
							vector<vector<int>> RemainedDimensionValue;
							RemainedDimensionValue.insert(RemainedDimensionValue.begin(), tempRemainedDimensionValue.begin(), tempRemainedDimensionValue.end());
							RemainedDimensionValue[k].insert(RemainedDimensionValue[k].begin(), expandDimensionValue.begin(), expandDimensionValue.end());
							mutantsDimensionValue.push_back(RemainedDimensionValue);
							//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
							for (int num = 0; num<expandCell.size(); num++)
							{
								matrixRow.push_back(expandCell[num]);
							}
							for (int num = 0; num<removedCells.size(); num++)
							{
								matrixRow.push_back(removedCells[num]);
							}
							corresponding_tests.push_back(matrixRow);
						}

						break;
					case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2��
						//���ǵÿ���k==iʱ��expand֮����ܺ�shrink֮ǰ��һ��
						int num = tempRemainedDimensionValue.at(k).at(0);
						bool noPlus1 = false;
						bool noMinus1 = false;
						if (k == i)
						{

							if (biggestValue == 2)//���ʱ���ܽ���expand
							{
								noPlus1 = true;
								noMinus1 = true;
							}
							else
							{
								int anotherNum = removedDimensionValues.at(j).at(0);
								if (commonFuncs.boundValue(num + 1, biggestValue) == anotherNum)
									noPlus1 = true;
								if (commonFuncs.boundValue(num - 1, biggestValue) == anotherNum)
									noMinus1 = true;
							}
						}
						else
						{
							if (biggestValue == 2)//���ʱ��ֻ�ܵ�����expand
							{
								if (num == 0)//�����ά���ֵ��2����ônumֻ������0��1
									noMinus1 = true;
								else
									noPlus1 = true;
							}
						}
						if (!noPlus1)
						{
							expandDimensionValue.clear();
							tempExpandDimensionValue.clear();
							expandCell.clear();
							expandDimensionValue.push_back(commonFuncs.boundValue(num + 1, biggestValue));
							commonFuncs.getTempDimensionValue(tempRemainedDimensionValue, k, expandDimensionValue, tempExpandDimensionValue);
							commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
							commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
							if (removedCells.size() + expandCell.size()>0){
								expandCells.push_back(expandCell);
								literalNum.push_back(literalN);
							}
							literalN++;//

							vector<CCell*> matrixRow;
							// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
							vector<vector<int>> RemainedDimensionValue;
							RemainedDimensionValue.insert(RemainedDimensionValue.begin(), tempRemainedDimensionValue.begin(), tempRemainedDimensionValue.end());
							RemainedDimensionValue[k].insert(RemainedDimensionValue[k].begin(), expandDimensionValue.begin(), expandDimensionValue.end());
							mutantsDimensionValue.push_back(RemainedDimensionValue);
							//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
							//corresponding_tests.push_back(expandCell);
							for (int num = 0; num<expandCell.size(); num++)
							{
								matrixRow.push_back(expandCell[num]);
							}
							for (int num = 0; num<removedCells.size(); num++)
							{
								matrixRow.push_back(removedCells[num]);
							}
							corresponding_tests.push_back(matrixRow);

						}
						if (!noMinus1)
						{
							expandDimensionValue.clear();
							tempExpandDimensionValue.clear();
							expandCell.clear();
							expandDimensionValue.push_back(commonFuncs.boundValue(num - 1, biggestValue));
							commonFuncs.getTempDimensionValue(tempRemainedDimensionValue, k, expandDimensionValue, tempExpandDimensionValue);
							commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
							commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
							if (removedCells.size() + expandCell.size()>0){
								expandCells.push_back(expandCell);
								literalNum.push_back(literalN);
							}
							literalN++;//

							vector<CCell*> matrixRow;
							// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
							vector<vector<int>> RemainedDimensionValue;
							RemainedDimensionValue.insert(RemainedDimensionValue.begin(), tempRemainedDimensionValue.begin(), tempRemainedDimensionValue.end());
							RemainedDimensionValue[k].insert(RemainedDimensionValue[k].begin(), expandDimensionValue.begin(), expandDimensionValue.end());
							mutantsDimensionValue.push_back(RemainedDimensionValue);
							//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
							//corresponding_tests.push_back(expandCell);
							for (int num = 0; num<expandCell.size(); num++)
							{
								matrixRow.push_back(expandCell[num]);
							}
							for (int num = 0; num<removedCells.size(); num++)
							{
								matrixRow.push_back(removedCells[num]);
							}
							corresponding_tests.push_back(matrixRow);
						}
						break;
					}

				}
				//�Ѿ���ȡ���е�expand cells����ȥ����һ��Ͷ�Ӧ�������һ��ϲ���ɾ����е�һ��
				for (int m = 0; m<expandCells.size(); m++)
				{
					vector<CCell*> matrixRow;
					for (int a = 0; a<removedCells.size(); a++)
						matrixRow.push_back(removedCells.at(a));
					for (int b = 0; b<expandCells.at(m).size(); b++)
						matrixRow.push_back(expandCells.at(m).at(b));
					if (param != 1)
						matrixCells.push_back(matrixRow);
				}
			}
		}
	}

	return "";
}
//����LRF�㷨
string LRFGenerateCellsForMutants(
	CCommonFuncs commonFuncs,
	string dimensionValuesStr, 
	bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests,
	vector<int> allOTPs, 
	vector<int> allUTPs
	)
{
	vector<int>literalNum;
	return runWithParameter(commonFuncs,0, literalNum, dimensionValuesStr, vOdd, matrixCells,
	 mutantsDimensionValue,
	 corresponding_tests, allOTPs, allUTPs);
}

