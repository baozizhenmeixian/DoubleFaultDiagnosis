#include "stdafx.h"
#include "DTLIFhandler.h"

//���������double fault lifxlif ����
string DTLIFhandler::run2(string dimensionValuesStr1, string dimensionValuesStr2, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> all3OTPs)
{
	vector<vector<int>> dimensionValues1;
	vector<int> biggestDimensionValueAllowed1;
	//��ȡ��άֵ
	m_commonFuncs.getDimensionValue(dimensionValuesStr1,false,dimensionValues1,biggestDimensionValueAllowed1);

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
					vector<vector<int>> dimensionValues2;
					vector<int> biggestDimensionValueAllowed2;
					//��ȡ��άֵ
					m_commonFuncs.getDimensionValue(dimensionValuesStr2,false,dimensionValues2,biggestDimensionValueAllowed2);

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
									shrinkCells.push_back(shrinkCell);
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
									shrinkCells.push_back(shrinkCell);
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
					//��ʱΪ��Ч��������Ҫ���ǵڶ��������
					validshrink = false;
					vector<vector<int>> dimensionValues2;
					vector<int> biggestDimensionValueAllowed2;
					//��ȡ��άֵ
					m_commonFuncs.getDimensionValue(dimensionValuesStr2,false,dimensionValues2,biggestDimensionValueAllowed2);

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
									shrinkCells.push_back(shrinkCell);
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
									shrinkCells.push_back(shrinkCell);
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
	



	//��ʾϵ������
	//showStr = getShowStr(matrixCells);
	//����ڴ�
	/*for (int i = 0; i < shrinkCells.size(); i++)
	{
		m_commonFuncs.clearCellList(shrinkCells.at(i));

	}*/
	shrinkCells.clear();
	vector<vector<CCell*>>().swap(shrinkCells);

	return "";

}

