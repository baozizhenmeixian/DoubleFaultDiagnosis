#include "stdafx.h"
#include "DTDFCommonFuncs.h"


DTDFCommonFuncs::DTDFCommonFuncs(void)
{
}


DTDFCommonFuncs::~DTDFCommonFuncs(void)
{
}

//���������double fault lifxlif ����
string DTDFCommonFuncs::generateMatrixofOTP(string dimensionValuesStr1, string dimensionValuesStr2,vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> all3OTPs,bool &hasequal, bool vOdd)
{
	vector<vector<int>> dimensionValues1;
	vector<int> biggestDimensionValueAllowed1;
	//��ȡ��άֵ
	m_commonFuncs.getDimensionValue(dimensionValuesStr1,vOdd,dimensionValues1,biggestDimensionValueAllowed1);

	int dimensionNum1 = dimensionValues1.size();
	vector<vector<CCell*>> shrinkCells;
	//��ÿһά������С����ȡ��С���ֵĵ�
	for(int i=0; i<dimensionNum1; i++)
	{
		vector<CCell*> shrinkCell; //�����ĸ�
		vector<CCell*> shrinkOtpCell1;//��һ��������otp
		vector<vector<int>> tempShrinkDimensionValue;
		vector<int> shrinkDimensionValue1;//��һ����ά���ϵ�����
		bool validshrink = false; //�ж���һά���ϵ������Ƿ���Ч
		int ithDimensionNum1 = dimensionValues1.at(i).size();
		switch(ithDimensionNum1)
		{
		case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
			for(int m=1; m<5; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkOtpCell1.clear();
				shrinkDimensionValue1.clear();
				shrinkDimensionValue1.push_back(m_commonFuncs.boundValue(m,4));
				shrinkDimensionValue1.push_back(m_commonFuncs.boundValue(m+1,4));
				m_commonFuncs.getTempDimensionValue(dimensionValues1,i,shrinkDimensionValue1,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
				m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell1,allOTPs);
				m_commonFuncs.kickOutCells(shrinkCell,allOTPs);
				if(shrinkCell.size()>0)
				{
					//shrinkCells.push_back(shrinkCell);
					validshrink = true;
					//m_commonFuncs.clearCellList(shrinkCell);
					m_commonFuncs.clearCellList(shrinkOtpCell1);

				}
				else
				{
					//��ʱΪ��Ч��������Ҫ���ǵڶ��������
					validshrink = false;
					hasequal = true;
					vector<vector<int>> dimensionValues2;
					vector<int> biggestDimensionValueAllowed2;
					//��ȡ��άֵ
					m_commonFuncs.getDimensionValue(dimensionValuesStr2,vOdd,dimensionValues2,biggestDimensionValueAllowed2);

					int dimensionNum2 = dimensionValues2.size();
					//��ʼö�ٵڶ������ά���ϵ�����
					for (int j = 0; j <dimensionNum2; j++)
					{
						vector<CCell*> shrinkOtpCell2;//�ڶ���������otp
						vector<CCell*> OtpCellsof2Terms;//���ͬ������otp
						vector<int> shrinkDimensionValue2;//�ڶ�����ά���ϵ�����
						int ithDimensionNum2 = dimensionValues2.at(j).size();
						switch(ithDimensionNum2)
						{
						case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
							for(int m=1; m<5; m++)
							{
								tempShrinkDimensionValue.clear();
								shrinkCell.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
								m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell,allOTPs);
								if(shrinkCell.size()>0)
								{
									//shrinkCells.push_back(shrinkCell);
									//m_commonFuncs.clearCellList(shrinkCell);
									m_commonFuncs.clearCellList(shrinkOtpCell2);
								}
								else
								{
									OtpCellsof2Terms.clear();
									m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
									m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);
									if (!OtpCellsof2Terms.empty())
									{
										shrinkCells.push_back(OtpCellsof2Terms);
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										//m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
									else
									{
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
								}
							}
							break;
						case 2: //�����������Ϊ1
							for(int m=0; m<2; m++)
							{
								tempShrinkDimensionValue.clear();
								shrinkCell.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(j).at(m),4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
								m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell,allOTPs);//Shrink�ĵ㲻����OTP
								if(shrinkCell.size()>0)
								{
									//shrinkCells.push_back(shrinkCell);
									//m_commonFuncs.clearCellList(shrinkCell);
									m_commonFuncs.clearCellList(shrinkOtpCell2);
								}
								else
								{
									OtpCellsof2Terms.clear();
									m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
									m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);
									if (!OtpCellsof2Terms.empty())
									{
										shrinkCells.push_back(OtpCellsof2Terms);
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										//m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
									else
									{
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
								}
							}
							break;
						case 1://������������ٽ�������
							break;
						}			
					}//�ڶ�����������

				}

			}
			break;
		case 2: //�����������Ϊ1
			for(int m=0; m<2; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkOtpCell1.clear();
				shrinkDimensionValue1.clear();
				shrinkDimensionValue1.push_back(m_commonFuncs.boundValue(dimensionValues1.at(i).at(m),4));
				m_commonFuncs.getTempDimensionValue(dimensionValues1,i,shrinkDimensionValue1,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
				m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell1,allOTPs);
				m_commonFuncs.kickOutCells(shrinkCell,allOTPs);//Shrink�ĵ㲻����OTP
				if(shrinkCell.size()>0)
				{
					//shrinkCells.push_back(shrinkCell);
					validshrink = true;
					//m_commonFuncs.clearCellList(shrinkCell);
					m_commonFuncs.clearCellList(shrinkOtpCell1);
				}
				else
				{
					validshrink = false;
					hasequal = true;
					//��ʱΪ��Ч��������Ҫ���ǵڶ��������
					vector<vector<int>> dimensionValues2;
					vector<int> biggestDimensionValueAllowed2;
					//��ȡ��άֵ
					m_commonFuncs.getDimensionValue(dimensionValuesStr2,vOdd,dimensionValues2,biggestDimensionValueAllowed2);

					int dimensionNum2 = dimensionValues2.size();
					//��ʼö�ٵڶ������ά���ϵ�����
					for (int j = 0; j <dimensionNum2; j++)
					{
						vector<CCell*> shrinkOtpCell2;//�ڶ���������otp
						vector<CCell*> OtpCellsof2Terms;//���ͬ������otp
						vector<int> shrinkDimensionValue2;//�ڶ�����ά���ϵ�����
						int ithDimensionNum2 = dimensionValues2.at(j).size();
						switch(ithDimensionNum2)
						{
						case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
							for(int m=1; m<5; m++)
							{
								tempShrinkDimensionValue.clear();
								shrinkCell.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
								m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell,allOTPs);
								if(shrinkCell.size()>0)
								{
									//shrinkCells.push_back(shrinkCell);
									//m_commonFuncs.clearCellList(shrinkCell);
									m_commonFuncs.clearCellList(shrinkOtpCell2);
								}
								else
								{
									OtpCellsof2Terms.clear();
									m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
									m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);
									if (!OtpCellsof2Terms.empty())
									{
										shrinkCells.push_back(OtpCellsof2Terms);
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										//m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
									else
									{
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
								}
							}
							break;
						case 2: //�����������Ϊ1
							for(int m=0; m<2; m++)
							{
								tempShrinkDimensionValue.clear();
								vector<vector<int>>().swap(tempShrinkDimensionValue);
								shrinkCell.clear();
								vector<CCell*>().swap(shrinkCell);
								shrinkOtpCell2.clear();
								vector<CCell*>().swap(shrinkOtpCell2);
								shrinkDimensionValue2.clear();
								vector<int>().swap(shrinkDimensionValue2);
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(j).at(m),4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
								m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell,allOTPs);//Shrink�ĵ㲻����OTP
								if(shrinkCell.size()>0)
								{
									//shrinkCells.push_back(shrinkCell);
									//m_commonFuncs.clearCellList(shrinkCell);
									m_commonFuncs.clearCellList(shrinkOtpCell2);
								}
								else
								{
									OtpCellsof2Terms.clear();
									m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
									m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);
									if (!OtpCellsof2Terms.empty())
									{
										shrinkCells.push_back(OtpCellsof2Terms);
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										//m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
									else
									{
										m_commonFuncs.clearCellList(shrinkCell);
										m_commonFuncs.clearCellList(shrinkOtpCell2);
										m_commonFuncs.clearCellList(OtpCellsof2Terms);
									}
								}
							}
							break;
						case 1://������������ٽ�������
							break;
						}
						//����ڴ�
					}//�ڶ�����������
				}
				
			}
			break;
		case 1://������������ٽ�������
			break;
		}
		//����ڴ�
	}//��һ��ö�ٽ���

	//�Ѿ���ȡ���е�shrink cells�������п��ܵ�shrink���һ������һ�б�ʾһ��shrink
	for(int i=0; i<shrinkCells.size(); i++)
	{
		vector<CCell*> matrixRow;
		for(int j=0; j<shrinkCells.at(i).size(); j++)
			matrixRow.push_back(shrinkCells.at(i).at(j));
		matrixCells.push_back(matrixRow);
	}
	shrinkCells.clear();
	vector<vector<CCell*>>().swap(shrinkCells);

	return "";

}

//���������double fault lofxlof
string DTDFCommonFuncs::generateMatrixofRFP(string dimensionValuesStr1, vector<vector<CCell*>>& matrixCells, vector<int> allUTPs, vector<int> allOTPs, bool &hasequal,bool vOdd)
{
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;

	m_commonFuncs.getDimensionValue(dimensionValuesStr1, vOdd, dimensionValues, biggestDimensionValueAllowed);
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
					expandDimensionValue.push_back(m_commonFuncs.boundValue(minNum + 2, biggestValue));
					expandDimensionValue.push_back(m_commonFuncs.boundValue(minNum + 3, biggestValue));
				}
				m_commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				m_commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				m_commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if (expandCell.size()>0){//��Ч���ţ����ü�¼��
					/*printf("[%d]expandCell\n", __LINE__);
					for (int j = 0; j < expandCell.size(); j++){
						printf("[%d]", expandCell[j]->getTModelIndex());
					}
					printf("\n");*/
					expandCell.clear();
				}
				{//����������
					hasequal = true;
					// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
					vector<vector<int>> RemainedDimensionValues;
					RemainedDimensionValues.insert(RemainedDimensionValues.begin(), dimensionValues.begin(), dimensionValues.end());
					RemainedDimensionValues[i].insert(RemainedDimensionValues[i].begin(), expandDimensionValue.begin(), expandDimensionValue.end());

					int dimensionNum2 = RemainedDimensionValues.size();
					//��ʼö�ٱ仯�����ά���ϵ�����
					for (int j = 0; j <dimensionNum2; j++)
					{
						vector<CCell*> expandCell2;
						vector<vector<int>> tempExpandDimensionValue2;
						vector<int> expandDimensionValue2;
						int biggestValue2 = 0;
						if (vOdd && j == dimensionNum2 - 1)
							biggestValue2 = 2;
						else
							biggestValue2 = 4;
						int ithDimensionNum2 = RemainedDimensionValues.at(j).size();
						switch (ithDimensionNum2)
						{
						case 4: //����������ܽ�������
							break;
						case 2: //�����������Ϊ4
							if (biggestValue2 == 4)//�������2��������
							{
								int num1 = RemainedDimensionValues.at(j).at(0);
								int num2 = RemainedDimensionValues.at(j).at(1);
								if ((num1 == 1 && num2 == 4) || (num1 == 4 && num2 == 1))
								{
									expandDimensionValue2.push_back(2);
									expandDimensionValue2.push_back(3);
								}
								else
								{
									int minNum = num1<num2 ? num1 : num2;
									expandDimensionValue2.push_back(m_commonFuncs.boundValue(minNum + 2, biggestValue2));
									expandDimensionValue2.push_back(m_commonFuncs.boundValue(minNum + 3, biggestValue2));
								}
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, i, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size()>0){
									expandCells.push_back(expandCell2);//��¼����
								}
							}
							break;
						case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2
							int num = RemainedDimensionValues.at(j).at(0);
							bool noPlus1 = false, noMinus1 = false;
							if (biggestValue2 == 2)//�������ֻ����һ�����䣬num��ֵ��1��2
							{
								if (num == 1)
									noMinus1 = true;
								else
									noPlus1 = true;
							}

							if (!noPlus1)
							{
								expandDimensionValue2.push_back(m_commonFuncs.boundValue(num + 1, biggestValue2));
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, i, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size() > 0)
								{
									expandCells.push_back(expandCell2);
								}
							}

							if (!noMinus1)
							{
								expandDimensionValue2.clear();
								tempExpandDimensionValue2.clear();
								expandCell2.clear();
								expandDimensionValue2.push_back(m_commonFuncs.boundValue(num - 1, biggestValue2));
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, i, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size() > 0)
								{
									expandCells.push_back(expandCell2);
								}
							}
							break;
						}
						expandCell2.clear();
						tempExpandDimensionValue2.clear();
						expandDimensionValue2.clear();
					}//�ڶ������Ž���
				}
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
				expandDimensionValue.push_back(m_commonFuncs.boundValue(num + 1, biggestValue));
				m_commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				m_commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				m_commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if (expandCell.size() > 0)
				{
					/*printf("[%d]expandCell\n", __LINE__);
					for (int j = 0; j < expandCell.size(); j++){
						printf("[%d]", expandCell[j]->getTModelIndex());
					}
					printf("\n");*/
					expandCell.clear();
				}
				{//����������
					hasequal = true;
					// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
					vector<vector<int>> RemainedDimensionValues;
					RemainedDimensionValues.insert(RemainedDimensionValues.begin(), dimensionValues.begin(), dimensionValues.end());
					RemainedDimensionValues[i].insert(RemainedDimensionValues[i].begin(), expandDimensionValue.begin(), expandDimensionValue.end());

					int dimensionNum2 = RemainedDimensionValues.size();
					//��ʼö�ٱ仯�����ά���ϵ�����
					for (int j = 0; j <dimensionNum2; j++)
					{
						vector<CCell*> expandCell2;
						vector<vector<int>> tempExpandDimensionValue2;
						vector<int> expandDimensionValue2;
						int biggestValue2 = 0;
						if (vOdd && j == dimensionNum2 - 1)
							biggestValue2 = 2;
						else
							biggestValue2 = 4;
						int ithDimensionNum2 = RemainedDimensionValues.at(j).size();
						switch (ithDimensionNum2)
						{
						case 4: //����������ܽ�������
							break;
						case 2: //�����������Ϊ4
							if (biggestValue2 == 4)//�������2��������
							{
								int num1 = RemainedDimensionValues.at(j).at(0);
								int num2 = RemainedDimensionValues.at(j).at(1);
								if ((num1 == 1 && num2 == 4) || (num1 == 4 && num2 == 1))
								{
									expandDimensionValue2.push_back(2);
									expandDimensionValue2.push_back(3);
								}
								else
								{
									int minNum = num1<num2 ? num1 : num2;
									expandDimensionValue2.push_back(m_commonFuncs.boundValue(minNum + 2, biggestValue2));
									expandDimensionValue2.push_back(m_commonFuncs.boundValue(minNum + 3, biggestValue2));
								}
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, j, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size()>0)
									expandCells.push_back(expandCell2);//��¼����
							}
							break;
						case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2
							int num = RemainedDimensionValues.at(j).at(0);
							bool noPlus1 = false, noMinus1 = false;
							if (biggestValue2 == 2)//�������ֻ����һ�����䣬num��ֵ��1��2
							{
								if (num == 1)
									noMinus1 = true;
								else
									noPlus1 = true;
							}

							if (!noPlus1)
							{
								expandDimensionValue2.push_back(m_commonFuncs.boundValue(num + 1, biggestValue2));
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, j, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size() > 0)
								{
									expandCells.push_back(expandCell2);
								}
							}

							if (!noMinus1)
							{
								expandDimensionValue2.clear();
								tempExpandDimensionValue2.clear();
								expandCell2.clear();
								expandDimensionValue2.push_back(m_commonFuncs.boundValue(num - 1, biggestValue2));
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, j, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size() > 0)
								{
									expandCells.push_back(expandCell2);
								}
							}
							break;
						}
						expandCell2.clear();
						tempExpandDimensionValue2.clear();
						expandDimensionValue2.clear();
					}//�ڶ������Ž���
				}
			}

			if (!noMinus1)
			{
				expandDimensionValue.clear();
				tempExpandDimensionValue.clear();
				expandCell.clear();
				expandDimensionValue.push_back(m_commonFuncs.boundValue(num - 1, biggestValue));
				m_commonFuncs.getTempDimensionValue(dimensionValues, i, expandDimensionValue, tempExpandDimensionValue);
				m_commonFuncs.ListAllCells(0, tempExpandDimensionValue, expandCell);
				m_commonFuncs.kickOutCells(expandCell, allOTPs, allUTPs);//expand�ĵ㲻����TP
				if (expandCell.size() > 0)
				{
					/*printf("[%d]expandCell\n", __LINE__);
					for (int j = 0; j < expandCell.size();j++){
						printf("[%d]", expandCell[j]->getTModelIndex());
					}
					printf("\n");*/
					expandCell.clear();
				}
				{//����������
					hasequal = true;
					// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)							
					vector<vector<int>> RemainedDimensionValues;
					RemainedDimensionValues.insert(RemainedDimensionValues.begin(), dimensionValues.begin(), dimensionValues.end());
					RemainedDimensionValues[i].insert(RemainedDimensionValues[i].begin(), expandDimensionValue.begin(), expandDimensionValue.end());
					int dimensionNum2 = RemainedDimensionValues.size();
					//��ʼö�ٱ仯�����ά���ϵ�����
					for (int j = 0; j <dimensionNum2; j++)
					{
						vector<CCell*> expandCell2;
						vector<vector<int>> tempExpandDimensionValue2;
						vector<int> expandDimensionValue2;
						int biggestValue2 = 0;
						if (vOdd && j == dimensionNum2 - 1)
							biggestValue2 = 2;
						else
							biggestValue2 = 4;
						int ithDimensionNum2 = RemainedDimensionValues.at(j).size();
						switch (ithDimensionNum2)
						{
						case 4: //����������ܽ�������
							break;
						case 2: //�����������Ϊ4
							if (biggestValue2 == 4)//�������2��������
							{
								int num1 = RemainedDimensionValues.at(j).at(0);
								int num2 = RemainedDimensionValues.at(j).at(1);
								if ((num1 == 1 && num2 == 4) || (num1 == 4 && num2 == 1))
								{
									expandDimensionValue2.push_back(2);
									expandDimensionValue2.push_back(3);
								}
								else
								{
									int minNum = num1<num2 ? num1 : num2;
									expandDimensionValue2.push_back(m_commonFuncs.boundValue(minNum + 2, biggestValue2));
									expandDimensionValue2.push_back(m_commonFuncs.boundValue(minNum + 3, biggestValue2));
								}
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, j, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size()>0){
									expandCells.push_back(expandCell2);//��¼����
								}
							}
							break;
						case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2
							int num = RemainedDimensionValues.at(j).at(0);
							bool noPlus1 = false, noMinus1 = false;
							if (biggestValue2 == 2)//�������ֻ����һ�����䣬num��ֵ��1��2
							{
								if (num == 1)
									noMinus1 = true;
								else
									noPlus1 = true;
							}

							if (!noPlus1)
							{
								expandDimensionValue2.push_back(m_commonFuncs.boundValue(num + 1, biggestValue2));
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, j, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size() > 0)
								{
									expandCells.push_back(expandCell2);
								}
							}

							if (!noMinus1)
							{
								expandDimensionValue2.clear();
								tempExpandDimensionValue2.clear();
								expandCell2.clear();
								expandDimensionValue2.push_back(m_commonFuncs.boundValue(num - 1, biggestValue2));
								m_commonFuncs.getTempDimensionValue(RemainedDimensionValues, j, expandDimensionValue2, tempExpandDimensionValue2);
								m_commonFuncs.ListAllCells(0, tempExpandDimensionValue2, expandCell2);
								m_commonFuncs.kickOutCells(expandCell2, allOTPs, allUTPs);//expand�ĵ㲻����TP
								if (expandCell2.size() > 0)
								{
									expandCells.push_back(expandCell2);
								}
							}
							break;
						}
						expandCell2.clear();
						tempExpandDimensionValue2.clear();
						expandDimensionValue2.clear();
					}//�ڶ������Ž���
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
	expandCells.clear();
	vector<vector<CCell*>>().swap(expandCells);
	return "";

}
