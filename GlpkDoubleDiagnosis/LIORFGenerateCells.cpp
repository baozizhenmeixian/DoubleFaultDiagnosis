#pragma once
#include "stdafx.h"
#include "LIORFGenerateCells.h"
#include "time.h"
using namespace std;
//����TOF�㷨
string TOFGenerateCells(CCommonFuncs commonFuncs, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs)
{
	string showStr = "";
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;

	commonFuncs.getDimensionValue(dimensionValuesStr, vOdd, dimensionValues, biggestDimensionValueAllowed);

	vector<vector<int>> tempDimensionValues = dimensionValues;
	int dimensionNum = dimensionValues.size();
	vector<CCell*> TPCell;
	//ȡ����UTP
	for (int i = 0; i<dimensionNum; i++)
	{
		vector<int> TPDimensionValue;
		vector<vector<int>> tempTPDimensionValue;
		int ithDimensionNum = dimensionValues.at(i).size();
		for (int j = 0; j<ithDimensionNum; j++)
		{
			TPDimensionValue.push_back(dimensionValues.at(i).at(j));
		}
		commonFuncs.getTempDimensionValue(tempDimensionValues, i, TPDimensionValue, tempTPDimensionValue);
		tempDimensionValues = tempTPDimensionValue;
	}
	commonFuncs.ListAllCells(0, tempDimensionValues, TPCell);
	commonFuncs.kickOutCells(TPCell, allOTPs);

	matrixCells.push_back(TPCell);

	//��ʾϵ������
	//showStr = getShowStr(matrixCells);

	return showStr;
}

//����LOF�㷨
string LOFGenerateCells(CCommonFuncs commonFuncs, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allUTPs)
{
	string showStr = "";
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;

	commonFuncs.getDimensionValue(dimensionValuesStr, vOdd, dimensionValues, biggestDimensionValueAllowed);
	int totalLOFliteralsNum = 0;
	int dimensionNum = dimensionValues.size();
	vector<vector<CCell*>> expandCells;
	//��ÿһά���зŴ�һ������ȡ�Ŵ󲿷ֵĵ�
	for (int i = 0; i<dimensionNum; i++)
	{
		vector<CCell*> expandCell;
		vector<vector<int>> tempExpandDimensionValue;
		vector<int> expandDimensionValue;
		int biggestValue = biggestDimensionValueAllowed.at(i);//��ά��������ֵ��2��4��
		int ithDimensionNum = dimensionValues.at(i).size();
		switch (ithDimensionNum)
		{
		case 4: //����������ܽ�������
			break;
		case 2: //�����������Ϊ4
			if (biggestValue == 4)//�������2��������
			{
				totalLOFliteralsNum += 1;
				int num1 = dimensionValues.at(i).at(0);
				int num2 = dimensionValues.at(i).at(1);
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
				commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if (expandCell.size()>0)
					expandCells.push_back(expandCell);
			}
			break;
		case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2
			int num = dimensionValues.at(i).at(0);
			bool noPlus1 = false, noMinus1 = false;
			if (biggestValue == 2)//�������ֻ����һ�����䣬num��ֵ��1��2
			{
				if (num == 1)
					noMinus1 = true;
				else
					noPlus1 = true;
			}

			if (!noPlus1)
			{
				expandDimensionValue.push_back(commonFuncs.boundValue(num + 1, biggestValue));
				commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if (expandCell.size() > 0)
				{
					totalLOFliteralsNum += 1;
					expandCells.push_back(expandCell);
				}
			}

			if (!noMinus1)
			{
				expandDimensionValue.clear();
				tempExpandDimensionValue.clear();
				expandCell.clear();
				expandDimensionValue.push_back(commonFuncs.boundValue(num - 1, biggestValue));
				commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if (expandCell.size() > 0)
				{
					expandCells.push_back(expandCell);
					totalLOFliteralsNum += 1;
				}
			}
			break;
		}
	}

	//�Ѿ���ȡ���е�expand cells�������п��ܵ��������һ������һ�б�ʾһ������
	for (int i = 0; i<expandCells.size(); i++)
	{
		vector<CCell*> matrixRow;
		for (int j = 0; j<expandCells.at(i).size(); j++)
			matrixRow.push_back(expandCells.at(i).at(j));
		matrixCells.push_back(matrixRow);
	}

	//��ʾ��Ч�ܵ�Expand��
	string temp;
	temp.append("%d", totalLOFliteralsNum);
	showStr.append(temp);
	return showStr;
}

//����LRF�㷨  paramΪ1��ֻ����validShrink����  paramΪ2��ֻ����expand���� ����ֵ����������
string LRFrunWithParameter(CCommonFuncs commonFuncs, int param, vector<int>&literalNum, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allUTPs)
{
	string showStr = "";
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;
	int allShrinksCount = 0;
	int validShrinksCount = 0;


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
			/*
			*
			*ȥ��OTP֮������жϣ�
			*1)���removedCells�л��е㣬��˵����shrink��һ����Чshrink������������ĵ㡣
			*2)���removedCellsΪ�ռ�������Ҫ����������ĵ�
			*/

			//��ʣ���
			if (removedCells.size() > 0)
			{
				//һ��
				vector<CCell*> matrixRow;
				for (int m = 0; m<removedCells.size(); m++)
				{
					matrixRow.push_back(removedCells.at(m));
				}
				if (param != 2)
					matrixCells.push_back(matrixRow);
				validShrinksCount++;
				continue;
			}

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
						literalN++;//
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

	//��ʾ��ЧShrink�����ܵ�Shrink��
	string temp;
	temp.append("%d", allShrinksCount);
	showStr.append(temp + "\t");
	temp.append("%d", validShrinksCount);
	showStr.append(temp);
	return showStr;
}

string LRFGenerateCells(
	CCommonFuncs commonFuncs,
	string dimensionValuesStr,
	bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<int> allOTPs,
	vector<int> allUTPs
	)
{
	vector<int>literalNum;
	return LRFrunWithParameter(commonFuncs, 0, literalNum, dimensionValuesStr, vOdd, matrixCells, allOTPs, allUTPs);
}


string LIORFrunWithParameter(CCommonFuncs commonFuncs, bool separateNeed, vector<int>&literalNum, vector<vector<CCell*>>& lifCells, vector<vector<CCell*>>& lofCells, vector<vector<CCell*>>& lrfCells, string dimensionValuesStr, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allUTPs, bool isSingleLitral)
{
	string showStr = "";
	vector<vector<int>> dimensionValues;//dimensionValuesStr�����ݽṹ�� ��һ��vector��ά�� �ڶ����Ǹ�ά���ϵ�ֵ
	vector<int> biggestDimensionValueAllowed;
	bool needLOF = true; //ֻҪ��һ��invalid shrink��Ӧ����չ���Լ�����е�LOF ����Ҫ��������

	int totalLOFliteralsNum = 0;
	int maxValidLOFliteralsNum = 0;

	commonFuncs.getDimensionValue(dimensionValuesStr, vOdd, dimensionValues, biggestDimensionValueAllowed);

	bool needTOF = true;
	int dimensionNum = dimensionValues.size();
	//���ÿһά����Ϊԭ����1/2����ѡһά����Ϊԭ��������

	for (int i = 0; i<dimensionNum; i++)
	{
		//�ֱ��ȡʣ�µ�һ���ȥ����һ��
		vector<vector<int>> remainedDimensionValues;//��һ��vector����Ŀ����shrink����Ŀ���ڶ���vector���ÿһ��shrink�漰��ֵ(ά��ȷ��[1,4])
		vector<vector<int>> removedDimensionValues;//ͬ��
		vector<int> remainedTemp;
		vector<int> removedTemp;

		int ithDimensionValueNum = dimensionValues.at(i).size();//�ڸ�ά����ֵ�Ŀ��
		//�о���һά���ܳ��ֵ������������¼���� 
		switch (ithDimensionValueNum)
		{
		case 4://��һά��1,2,3,4�����
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

		//���ÿһ��shrink����**����**������Ϊ��������ȡԼ����ȥ����һ��������ĵ㣩
		for (int j = 0; j<removedDimensionValues.size(); j++)
		{
			//��ȡȥ������һ��ĵ㣬�����һ��ȥ���ж�������
			vector<CCell*> removedCells;
			vector<vector<int>> tempDimensionValue;

			//��ȡ��ȥ����һ��ĵ�,�������ĸ�ά��ֵ��Ŀǰ���ڵ�ά�ȣ�����������ֵ�����tempDimensionValue��������Ҳ�����
			commonFuncs.getTempDimensionValue(dimensionValues, i, removedDimensionValues.at(j), tempDimensionValue);
			commonFuncs.ListAllCells(0, tempDimensionValue, removedCells);//tempDimensionValue��һ����ά����ʽ��ֵ
			commonFuncs.kickOutCells(removedCells, allOTPs);//ȥ����һ���в���ѡ��OTP

			/*
			*
			*ȥ��OTP֮������жϣ�
			*1)���removedCells�л��е㣬��˵����shrink��һ����Чshrink������������ĵ㡣
			*2)���removedCellsΪ�ռ�������Ҫ����������ĵ�
			*/

			//��ʣ���:��һ�����
			if (removedCells.size() > 0)
			{
				//һ��
				vector<CCell*> matrixRow;
				for (int m = 0; m<removedCells.size(); m++)
				{
					matrixRow.push_back(removedCells.at(m));
				}
				matrixCells.push_back(matrixRow);
				needTOF = false;
				if (separateNeed)
					lifCells.push_back(matrixRow);//lif
				continue;
			}
			//needLOF = false;
			//��ʣ���:�ڶ������
			//��ȡʣ�µ�һ���dimensionֵ
			vector<vector<int>> tempRemainedDimensionValue;
			commonFuncs.getTempDimensionValue(dimensionValues, i, remainedDimensionValues.at(j), tempRemainedDimensionValue);

			//�������ȥ�����г����п��ܵ�����
			vector<vector<CCell*>> expandCells;
			int literalN = 0;//
			//�������ά�ȵ����䶼�ҵ�����Ч�ĵ㣬��˵����shrink��Ӧ������㱣֤�˸����LOF����Ҫ�ٵ�������
			bool needLOFInThisShrink = false;
			int validLOFliteralsNumThisShrink = 0;
			int totalLOFliteralsNumThisShrink = 0;
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
						totalLOFliteralsNumThisShrink += 1;
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
						if (expandCell.size()>0)
						{
							validLOFliteralsNumThisShrink += 1;
							expandCells.push_back(expandCell);
							literalNum.push_back(literalN);//
						}
						else
							needLOFInThisShrink = true;//�õ���literal��LOF����ͨ����shrink�⵽(��Ϊ2���Ա����ǵ���)
						literalN++;
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
						totalLOFliteralsNumThisShrink += 1;
						expandDimensionValue.clear();
						tempExpandDimensionValue.clear();
						expandCell.clear();
						expandDimensionValue.push_back(commonFuncs.boundValue(num + 1, biggestValue));
						commonFuncs.getTempDimensionValue(tempRemainedDimensionValue, k, expandDimensionValue, tempExpandDimensionValue);
						commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
						commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
						if (expandCell.size()>0)
						{
							validLOFliteralsNumThisShrink += 1;
							expandCells.push_back(expandCell);
							literalNum.push_back(literalN);//
						}
						else
							needLOFInThisShrink = true;//������literal�е�һ��literal��LOF����ͨ����shrink��⵽
						literalN++;
					}
					if (!noMinus1)
					{
						totalLOFliteralsNumThisShrink += 1;
						expandDimensionValue.clear();
						tempExpandDimensionValue.clear();
						expandCell.clear();
						expandDimensionValue.push_back(commonFuncs.boundValue(num - 1, biggestValue));
						commonFuncs.getTempDimensionValue(tempRemainedDimensionValue, k, expandDimensionValue, tempExpandDimensionValue);
						commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
						commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
						if (expandCell.size()>0)
						{
							validLOFliteralsNumThisShrink += 1;
							expandCells.push_back(expandCell);
							literalNum.push_back(literalN);
						}
						else
							needLOFInThisShrink = true;//������literal�е�һ��literal��LOF����ͨ����shrink��⵽
						literalN++;
					}
					break;
				}

			}

			//old:������һ��ά�ȵ�shrink�����invalid ����Ҫ����LOF ��Ϊexpand�ĵ�϶��ܼ�⵽LOF
			//new:������һ��ά�ȵ����䶼�ҵ��˿��õĵ㣬��˵����shrink���Լ������LOF
			if (!needLOFInThisShrink)
				needLOF = false;
			if (validLOFliteralsNumThisShrink > maxValidLOFliteralsNum)
				maxValidLOFliteralsNum = validLOFliteralsNumThisShrink;
			if (totalLOFliteralsNum == 0)
				totalLOFliteralsNum = totalLOFliteralsNumThisShrink;
			//�Ѿ���ȡ���е�expand cells����ȥ����һ��Ͷ�Ӧ�������һ��ϲ���ɾ����е�һ��(����Ľ�֮������Ӧ��ֻ����չ�ĵ�)
			for (int m = 0; m<expandCells.size(); m++)
			{
				vector<CCell*> matrixRow;
				for (int a = 0; a<removedCells.size(); a++)
					matrixRow.push_back(removedCells.at(a));
				for (int b = 0; b<expandCells.at(m).size(); b++)
					matrixRow.push_back(expandCells.at(m).at(b));
				matrixCells.push_back(matrixRow);
				if (separateNeed)
					lrfCells.push_back(matrixRow);//lrf
			}
		}
	}

	//��Ҫ��������LOF 
	if (!isSingleLitral && needLOF)
	{
		LOFGenerateCells(commonFuncs, dimensionValuesStr, vOdd, matrixCells, allOTPs, allUTPs);
		if (separateNeed)
			LOFGenerateCells(commonFuncs, dimensionValuesStr, vOdd, lofCells, allOTPs, allUTPs);//lof
	}

	if (needTOF)
	{//�����ֻ���һ���ת���TOF
		TOFGenerateCells(commonFuncs, dimensionValuesStr, vOdd, matrixCells, allOTPs);

	}
	if (isSingleLitral)
	{
		totalLOFliteralsNum = 0;
		maxValidLOFliteralsNum = 0;
	}
	//��ʾ��ЧExpand�����ܵ�Expand��
	string temp;
	temp.append("%d", totalLOFliteralsNum);
	showStr.append(temp + "\t");
	temp.append("%d", maxValidLOFliteralsNum);
	showStr.append(temp);
	return showStr;
}

string LIORFGenerateCells(CCommonFuncs commonFuncs, string dimensionValuesStr, bool vOdd, vector<int> allOTPs, vector<int> allUTPs, bool isSingleLitral, vector<vector<CCell*>>& matrixCells)
{
	vector<int>literalNum;
	vector<vector<CCell*>> lifCells;
	vector<vector<CCell*>> lofCells;
	vector<vector<CCell*>> lrfCells;
	string ret = LIORFrunWithParameter(commonFuncs, false, literalNum, lifCells, lofCells, lrfCells, dimensionValuesStr, vOdd, matrixCells, allOTPs, allUTPs, isSingleLitral);
	return ret;
}

string LIORFGenerateCells(CCommonFuncs commonFuncs, string dimensionValuesStr, bool vOdd, vector<int> allOTPs, vector<int> allUTPs, bool isSingleLitral, vector<vector<CCell*>>& matrixCells,
vector<vector<CCell*>> &lifCells,
vector<vector<CCell*>> &lofCells,
vector<vector<CCell*>> &lrfCells
)
{
	vector<int>literalNum;
	string ret = LIORFrunWithParameter(commonFuncs, true, literalNum, lifCells, lofCells, lrfCells, dimensionValuesStr, vOdd, matrixCells, allOTPs, allUTPs, isSingleLitral);
	return ret;
}