#include "stdafx.h"
#include "DTLIFhandler.h"


DTLIFhandler::DTLIFhandler(void)
{
}


DTLIFhandler::~DTLIFhandler(void)
{
}

// ����DTLIF�㷨����Լ������
string DTLIFhandler::run(string dimensionValuesStr, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs)
{
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;
	//��ȡ��άֵ
	m_commonFuncs.getDimensionValue(dimensionValuesStr,false,dimensionValues,biggestDimensionValueAllowed);

	int dimensionNum = dimensionValues.size();
	vector<vector<CCell*>> shrinkCells;
	//��ÿһά������С����ȡ��С���ֵĵ�
	for(int i=0; i<dimensionNum; i++)
	{
		vector<CCell*> shrinkCell;
		vector<vector<CCell*>> shrinkOtpCell;
		vector<vector<int>> tempShrinkDimensionValue;
		vector<int> shrinkDimensionValue;
		int ithDimensionNum = dimensionValues.at(i).size();
		switch(ithDimensionNum)
		{
		case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
			for(int m=1; m<5; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(m_commonFuncs.boundValue(m,4));
				shrinkDimensionValue.push_back(m_commonFuncs.boundValue(m+1,4));
				m_commonFuncs.getTempDimensionValue(dimensionValues,i,shrinkDimensionValue,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
				//m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell,allOTPs);
				m_commonFuncs.kickOutCells(shrinkCell,allOTPs);
				if(shrinkCell.size()>0)
					shrinkCells.push_back(shrinkCell);
				else
				{
					if(shrinkOtpCell.size()>0)
					{
						for(int i=0;i<shrinkOtpCell.size();i++)
						{
							shrinkCells.push_back(shrinkOtpCell.at(i));
						}
					}
				}
				
			}
			break;
		case 2: //�����������Ϊ1
			for(int m=0; m<2; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(m_commonFuncs.boundValue(dimensionValues.at(i).at(m),4));
				m_commonFuncs.getTempDimensionValue(dimensionValues,i,shrinkDimensionValue,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
				//m_commonFuncs.getShrinkOtps(shrinkCell,shrinkOtpCell,allOTPs);
				m_commonFuncs.kickOutCells(shrinkCell,allOTPs);//Shrink�ĵ㲻����OTP
				if(shrinkCell.size()>0)
					shrinkCells.push_back(shrinkCell);
				else
				{
					if(shrinkOtpCell.size()>0)
					{
						for(int i=0;i<shrinkOtpCell.size();i++)
						{
							shrinkCells.push_back(shrinkOtpCell.at(i));
						}
					}
				}
				
			}
			break;
		case 1://������������ٽ�������
			break;
		}
	}

	//�Ѿ���ȡ���е�shrink cells�������п��ܵ�shrink���һ������һ�б�ʾһ��shrink
	for(int i=0; i<shrinkCells.size(); i++)
	{
		vector<CCell*> matrixRow;
		for(int j=0; j<shrinkCells.at(i).size(); j++)
			matrixRow.push_back(shrinkCells.at(i).at(j));
		matrixCells.push_back(matrixRow);
	}

	//������opts׷�ӵ����Ͼ���һ�б�ʾһ��opt
	/*for(int i=0; i<allOTPs.size; i++)
	{
		vector<CCell*> matrixRow;
		matrixRow.push_back(all0OTPs.at(i));
		matrixCells.push_back(matrixRow);
	}*/

	//��ʾϵ������
	//showStr = getShowStr(matrixCells);

	return "";
}

string DTLIFhandler::run(
	string	dimensionValuesStr,
	vector<vector<CCell*>>& matrixCells,
	bool mutantsNeed, 
	vector<vector<vector<int>>>& mutantsDimensionValue,
	vector<vector<CCell*>>& corresponding_tests,
	vector<int> allOTPs
	)
{
	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;
	//��ȡ��άֵ
	m_commonFuncs.getDimensionValue(dimensionValuesStr,false,dimensionValues,biggestDimensionValueAllowed);

	int dimensionNum = dimensionValues.size();
	vector<vector<CCell*>> shrinkCells;
	//��ÿһά������С����ȡ��С���ֵĵ�
	for(int i=0; i<dimensionNum; i++)
	{
		vector<CCell*> shrinkCell;
		vector<vector<int>> tempShrinkDimensionValue;
		vector<int> shrinkDimensionValue;
		int ithDimensionNum = dimensionValues.at(i).size();
		switch(ithDimensionNum)
		{
		case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
			for(int m=1; m<5; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(m_commonFuncs.boundValue(m,4));
				shrinkDimensionValue.push_back(m_commonFuncs.boundValue(m+1,4));
				//tempShrinkDimensionValue  : �̺���ǰ�����ȱ�ݸ񼯺�
				m_commonFuncs.getTempDimensionValue(dimensionValues,i,shrinkDimensionValue,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
				m_commonFuncs.kickOutCells(shrinkCell,allOTPs);
				if(shrinkCell.size()>0){
					shrinkCells.push_back(shrinkCell);
				}
				if (mutantsNeed)
				{
					// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)
					shrinkDimensionValue.clear();
					shrinkDimensionValue.push_back(m_commonFuncs.boundValue(m+2,4));
					shrinkDimensionValue.push_back(m_commonFuncs.boundValue(m+3,4));
					vector<vector<int>> tempMutantDimensionValue;
					m_commonFuncs.getTempDimensionValue(dimensionValues,i,shrinkDimensionValue,tempMutantDimensionValue);
					mutantsDimensionValue.push_back(tempMutantDimensionValue);

					//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
					corresponding_tests.push_back(shrinkCell);
				}
			}
			break;
		case 2: //�����������Ϊ1
			for(int m=0; m<2; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkDimensionValue.clear();
				shrinkDimensionValue.push_back(m_commonFuncs.boundValue(dimensionValues.at(i).at(m),4));
				m_commonFuncs.getTempDimensionValue(dimensionValues,i,shrinkDimensionValue,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell);
				m_commonFuncs.kickOutCells(shrinkCell,allOTPs);//Shrink�ĵ㲻����OTP
				if(shrinkCell.size()>0){
					shrinkCells.push_back(shrinkCell);
				}
				if (mutantsNeed)
				{
					// ��ǰ��mutant,(��һ��������ʾ����������ά��ȡֵ����mutantsDimensionValue��)
					shrinkDimensionValue.clear();
					shrinkDimensionValue.push_back(m_commonFuncs.boundValue(dimensionValues.at(i).at((m+1)%2),4));
					vector<vector<int>> tempMutantDimensionValue;
					m_commonFuncs.getTempDimensionValue(dimensionValues,i,shrinkDimensionValue,tempMutantDimensionValue);
					mutantsDimensionValue.push_back(tempMutantDimensionValue);

					//�ѵ�ǰmutant��Ӧ�Ĳ����������ռ���corresponding_tests��
					corresponding_tests.push_back(shrinkCell);
				}
			}
			break;
		case 1://������������ٽ�������
			break;
		}
	}

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

	return "";
}
//���������double fault lifxlif ����
string DTLIFhandler::run(string dimensionValuesStr1, string dimensionValuesStr2, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> all3OTPs)
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
					shrinkCells.push_back(shrinkCell);
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
					shrinkCells.push_back(shrinkCell);
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
//����
string DTLIFhandler::run(string dimensionValuesStr1, string dimensionValuesStr2, vector<vector<CCell*>>& matrixCells, bool mutantsNeed, vector<vector<vector<int>>>& mutantsDimensionValue, vector<vector<CCell*>>& corresponding_tests, vector<int> allOTPs, vector<int> all3OTPs)
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
		vector<CCell*> shrinkCell1; //��һ�������ĸ�
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
				shrinkCell1.clear();
				shrinkOtpCell1.clear();
				shrinkDimensionValue1.clear();
				shrinkDimensionValue1.push_back(m_commonFuncs.boundValue(m,4));
				shrinkDimensionValue1.push_back(m_commonFuncs.boundValue(m+1,4));
				m_commonFuncs.getTempDimensionValue(dimensionValues1,i,shrinkDimensionValue1,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell1);
				m_commonFuncs.getShrinkOtps(shrinkCell1,shrinkOtpCell1,allOTPs);
				m_commonFuncs.kickOutCells(shrinkCell1,allOTPs);

				if(shrinkCell1.size()>0)
				{
					for(int u=0; u<shrinkCell1.size(); u++)
					{
						shrinkCell.push_back(shrinkCell1[u]);
					}
				}

				//���ǵڶ��������
				vector<vector<int>> dimensionValues2;
				vector<int> biggestDimensionValueAllowed2;
				//��ȡ��άֵ
				m_commonFuncs.getDimensionValue(dimensionValuesStr2,false,dimensionValues2,biggestDimensionValueAllowed2);

				int dimensionNum2 = dimensionValues2.size();
				//��ʼö�ٵڶ������ά���ϵ�����
				for (int j = 0; j <dimensionNum2; j++)
				{
					vector<CCell*> shrinkCell2; //�ڶ��������ĸ�
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
								shrinkCell2.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell2);
								m_commonFuncs.getShrinkOtps(shrinkCell2,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);

								if(shrinkCell2.size()>0)
								{
									for(int v=0; v<shrinkCell2.size(); v++)
									{
										shrinkCell.push_back(shrinkCell2[v]);
									}
								}

								OtpCellsof2Terms.clear();
								m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
								m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);

								if (!OtpCellsof2Terms.empty())
								{
									for(int w=0; w<OtpCellsof2Terms.size(); w++)
									{
										shrinkCell.push_back(OtpCellsof2Terms[w]);
									}
								}
								if(!shrinkCell.empty())
									shrinkCells.push_back(shrinkCell);
								/*m_commonFuncs.clearCellList(shrinkCell);
								m_commonFuncs.clearCellList(shrinkCell1);
								m_commonFuncs.clearCellList(shrinkCell2);
								m_commonFuncs.clearCellList(shrinkOtpCell1);
								m_commonFuncs.clearCellList(shrinkOtpCell2);
								m_commonFuncs.clearCellList(OtpCellsof2Terms);*/
							}
							break;
						case 2: //�����������Ϊ1
							for(int m=0; m<2; m++)
							{
								tempShrinkDimensionValue.clear();
								shrinkCell2.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(j).at(m),4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell2);
								m_commonFuncs.getShrinkOtps(shrinkCell2,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);

								if(shrinkCell2.size()>0)
								{
									for(int v=0; v<shrinkCell2.size(); v++)
									{
										shrinkCell.push_back(shrinkCell2[v]);
									}
								}

								OtpCellsof2Terms.clear();
								m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
								m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);

								if (!OtpCellsof2Terms.empty())
								{
									for(int w=0; w<OtpCellsof2Terms.size(); w++)
									{
										shrinkCell.push_back(OtpCellsof2Terms[w]);
									}
								}

								if(!shrinkCell.empty())
									shrinkCells.push_back(shrinkCell);
								/*m_commonFuncs.clearCellList(shrinkCell);
								m_commonFuncs.clearCellList(shrinkCell1);
								m_commonFuncs.clearCellList(shrinkCell2);
								m_commonFuncs.clearCellList(shrinkOtpCell1);
								m_commonFuncs.clearCellList(shrinkOtpCell2);
								m_commonFuncs.clearCellList(OtpCellsof2Terms);*/
							}
							break;
						case 1://������������ٽ�������
							break;			
					}
				}//�ڶ�����������
			}
			break;
		case 2://�����������Ϊ1
			for(int m=0; m<2; m++)
			{
				tempShrinkDimensionValue.clear();
				shrinkCell.clear();
				shrinkCell1.clear();
				shrinkOtpCell1.clear();
				shrinkDimensionValue1.clear();
				shrinkDimensionValue1.push_back(m_commonFuncs.boundValue(dimensionValues1.at(i).at(m),4));
				m_commonFuncs.getTempDimensionValue(dimensionValues1,i,shrinkDimensionValue1,tempShrinkDimensionValue);
				m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell1);
				m_commonFuncs.getShrinkOtps(shrinkCell1,shrinkOtpCell1,allOTPs);
				m_commonFuncs.kickOutCells(shrinkCell1,allOTPs);

				if(shrinkCell1.size()>0)
				{
					for(int u=0; u<shrinkCell1.size(); u++)
					{
						shrinkCell.push_back(shrinkCell1[u]);
					}
				}

				//���ǵڶ��������
				vector<vector<int>> dimensionValues2;
				vector<int> biggestDimensionValueAllowed2;
				//��ȡ��άֵ
				m_commonFuncs.getDimensionValue(dimensionValuesStr2,false,dimensionValues2,biggestDimensionValueAllowed2);

				int dimensionNum2 = dimensionValues2.size();
				//��ʼö�ٵڶ������ά���ϵ�����
				for (int j = 0; j <dimensionNum2; j++)
				{
					vector<CCell*> shrinkCell2; //�ڶ��������ĸ�
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
								shrinkCell2.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell2);
								m_commonFuncs.getShrinkOtps(shrinkCell2,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);

								if(shrinkCell2.size()>0)
								{
									for(int v=0; v<shrinkCell2.size(); v++)
									{
										shrinkCell.push_back(shrinkCell2[v]);
									}
								}

								OtpCellsof2Terms.clear();
								m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
								m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);

								if (!OtpCellsof2Terms.empty())
								{
									for(int w=0; w<OtpCellsof2Terms.size(); w++)
									{
										shrinkCell.push_back(OtpCellsof2Terms[w]);
									}
								}

								if(!shrinkCell.empty())
									shrinkCells.push_back(shrinkCell);
								/*m_commonFuncs.clearCellList(shrinkCell);
								m_commonFuncs.clearCellList(shrinkCell1);
								m_commonFuncs.clearCellList(shrinkCell2);
								m_commonFuncs.clearCellList(shrinkOtpCell1);
								m_commonFuncs.clearCellList(shrinkOtpCell2);
								m_commonFuncs.clearCellList(OtpCellsof2Terms);*/
							}
							break;
						case 2: //�����������Ϊ1
							for(int m=0; m<2; m++)
							{
								tempShrinkDimensionValue.clear();
								shrinkCell2.clear();
								shrinkOtpCell2.clear();
								shrinkDimensionValue2.clear();
								shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(j).at(m),4));
								m_commonFuncs.getTempDimensionValue(dimensionValues2,j,shrinkDimensionValue2,tempShrinkDimensionValue);
								m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue,shrinkCell2);
								m_commonFuncs.getShrinkOtps(shrinkCell2,shrinkOtpCell2,allOTPs);
								m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);

								if(shrinkCell2.size()>0)
								{
									for(int v=0; v<shrinkCell2.size(); v++)
									{
										shrinkCell.push_back(shrinkCell2[v]);
									}
								}

								OtpCellsof2Terms.clear();
								m_commonFuncs.getOverlapCells(shrinkOtpCell1,shrinkOtpCell2,OtpCellsof2Terms);
								m_commonFuncs.kickOutCells(OtpCellsof2Terms,all3OTPs);

								if (!OtpCellsof2Terms.empty())
								{
									for(int w=0; w<OtpCellsof2Terms.size(); w++)
									{
										shrinkCell.push_back(OtpCellsof2Terms[w]);
									}
								}

								if(!shrinkCell.empty())
									shrinkCells.push_back(shrinkCell);
								/*m_commonFuncs.clearCellList(shrinkCell);
								m_commonFuncs.clearCellList(shrinkCell1);
								m_commonFuncs.clearCellList(shrinkCell2);
								m_commonFuncs.clearCellList(shrinkOtpCell1);
								m_commonFuncs.clearCellList(shrinkOtpCell2);
								m_commonFuncs.clearCellList(OtpCellsof2Terms);*/
							}
							break;
						case 1://������������ٽ�������
							break;			
					}
				}//�ڶ�����������
			}
			break;
		case 1://������������ٽ�������
			break;
		}
	}//��һ��ö�ٽ���

	//�Ѿ���ȡ���е�shrink cells�������п��ܵ�shrink���һ������һ�б�ʾһ��shrink
	for(int i=0; i<shrinkCells.size(); i++)
	{
		matrixCells.push_back(shrinkCells[i]);
	}

	shrinkCells.clear();
	vector<vector<CCell*>>().swap(shrinkCells);

	return "";
}

