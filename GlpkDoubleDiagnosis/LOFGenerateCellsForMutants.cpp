#include "StdAfx.h"
#include "LOFGenerateCellsForMutants.h"

string LOFGenerateCellsForMutants(
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
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;

	commonFuncs.getDimensionValue(dimensionValuesStr, vOdd, dimensionValues, biggestDimensionValueAllowed);
	int totalLOFliteralsNum = 0;
	int dimensionNum = dimensionValues.size();
	vector<vector<CCell*>> expandCells;
	//��ÿһά���зŴ�һ������ȡ�Ŵ󲿷ֵĵ�
	for(int i=0; i<dimensionNum; i++)
	{
		vector<CCell*> expandCell;
		vector<vector<int>> tempExpandDimensionValue;
		vector<int> expandDimensionValue;
		int biggestValue = biggestDimensionValueAllowed.at(i);//��ά��������ֵ��2��4��
		int ithDimensionNum = dimensionValues.at(i).size();
		switch(ithDimensionNum)
		{
		case 4: //����������ܽ�������
			break;
		case 2: //�����������Ϊ4
			if(biggestValue==4)//�������2��������
			{
				totalLOFliteralsNum += 1;
				int num1 = dimensionValues.at(i).at(0);
				int num2 = dimensionValues.at(i).at(1);
				if((num1==1&&num2==4) || (num1==4&&num2==1))
				{
					expandDimensionValue.push_back(2);
					expandDimensionValue.push_back(3);
				}
				else
				{
					int minNum = num1<num2? num1:num2;
					expandDimensionValue.push_back(commonFuncs.boundValue(minNum + 2, biggestValue));
					expandDimensionValue.push_back(commonFuncs.boundValue(minNum + 3, biggestValue));
				}
				commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				/*printf("---------------tempShrinkDimensionValue start--------------\n");
				for (int z = 0; z<tempExpandDimensionValue.size(); z++){
					for (int y = 0; y < tempExpandDimensionValue[z].size(); y++){
						printf("%d ", tempExpandDimensionValue[z][y]);
					}
					printf("\n");
				}
				printf("\n");
				printf("---------------tempShrinkDimensionValue end--------------\n");*/
				commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if(expandCell.size()>0)
					expandCells.push_back(expandCell);

				// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
				vector<vector<int>> RemainedDimensionValue;
				RemainedDimensionValue.insert(RemainedDimensionValue.begin(),dimensionValues.begin(),dimensionValues.end());
				RemainedDimensionValue[i].insert(RemainedDimensionValue[i].begin(),expandDimensionValue.begin(),expandDimensionValue.end());
				mutantsDimensionValue.push_back(RemainedDimensionValue);
				//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
				corresponding_tests.push_back(expandCell);
			}
			break;
		case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2
			int num = dimensionValues.at(i).at(0);
			bool noPlus1=false, noMinus1=false;
			if(biggestValue==2)//�������ֻ����һ�����䣬num��ֵ��1��2
			{
				if(num==1)
					noMinus1 = true;
				else
					noPlus1 = true;
			}
			
			if(!noPlus1)
			{
				expandDimensionValue.push_back(commonFuncs.boundValue(num + 1, biggestValue));
				commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				/*printf("---------------tempShrinkDimensionValue start--------------\n");
				for (int z = 0; z<tempExpandDimensionValue.size(); z++){
					for (int y = 0; y < tempExpandDimensionValue[z].size(); y++){
						printf("%d ", tempExpandDimensionValue[z][y]);
					}
					printf("\n");
				}
				printf("\n");
				printf("---------------tempShrinkDimensionValue end--------------\n");*/
				commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if(expandCell.size() > 0 )
				{
					totalLOFliteralsNum += 1;
					expandCells.push_back(expandCell);
				}

				// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
				vector<vector<int>> RemainedDimensionValue;
				RemainedDimensionValue.insert(RemainedDimensionValue.begin(), dimensionValues.begin(), dimensionValues.end());
				RemainedDimensionValue[i].insert(RemainedDimensionValue[i].begin(), expandDimensionValue.begin(), expandDimensionValue.end());
				mutantsDimensionValue.push_back(RemainedDimensionValue);
				//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
				corresponding_tests.push_back(expandCell);
			}

			if(!noMinus1)
			{
				expandDimensionValue.clear();
				tempExpandDimensionValue.clear();
				expandCell.clear();
				expandDimensionValue.push_back(commonFuncs.boundValue(num - 1, biggestValue));
				commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				/*printf("---------------tempShrinkDimensionValue start--------------\n");
				for (int z = 0; z<tempExpandDimensionValue.size(); z++){
					for (int y = 0; y < tempExpandDimensionValue[z].size(); y++){
						printf("%d ", tempExpandDimensionValue[z][y]);
					}
					printf("\n");
				}
				printf("\n");
				printf("---------------tempShrinkDimensionValue end--------------\n");*/
				commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if(expandCell.size() > 0)
				{
					expandCells.push_back(expandCell);
					totalLOFliteralsNum += 1;
				}

				// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
				vector<vector<int>> RemainedDimensionValue;
				RemainedDimensionValue.insert(RemainedDimensionValue.begin(), dimensionValues.begin(), dimensionValues.end());
				RemainedDimensionValue[i].insert(RemainedDimensionValue[i].begin(), expandDimensionValue.begin(), expandDimensionValue.end());
				mutantsDimensionValue.push_back(RemainedDimensionValue);
				//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
				corresponding_tests.push_back(expandCell);
			}
			break;
		}
	}

	//�Ѿ���ȡ���е�expand cells�������п��ܵ��������һ������һ�б�ʾһ������
	for(int i=0; i<expandCells.size(); i++)
	{
		vector<CCell*> matrixRow;
		for(int j=0; j<expandCells.at(i).size(); j++)
			matrixRow.push_back(expandCells.at(i).at(j));
		matrixCells.push_back(matrixRow);
	}

	return "";
}
