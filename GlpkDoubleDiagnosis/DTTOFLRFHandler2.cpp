#include "stdafx.h"
#include "DTTOFLRFHandler.h"


//����DLRF�㷨  
string CDTTOFLRFHandler::run2(string dimensionValuesStr1, string dimensionValuesStr2, bool vOdd, vector<vector<CCell*>>& matrixCells, vector<int> allOTPs, vector<int> allUTPs)
{
	vector<vector<int>> dimensionValues1;
	vector<int> biggestDimensionValueAllowed1;


	m_commonFuncs.getDimensionValue(dimensionValuesStr1,vOdd,dimensionValues1,biggestDimensionValueAllowed1);

	int dimensionNum1 = dimensionValues1.size();
	//���ÿһά����Ϊԭ����1/2����ѡһά����Ϊԭ��������

	for(int i=0; i<dimensionNum1; i++)
	{
		//�ֱ��ȡʣ�µ�һ���ȥ����һ��
		vector<vector<int>> remainedDimensionValues1;
		vector<vector<int>> removedDimensionValues1;
		vector<int> remainedTemp1;
		vector<int> removedTemp1;

		int ithDimensionValueNum1 = dimensionValues1.at(i).size();
		//�о���һά���ܳ��ֵ������������¼����
		switch(ithDimensionValueNum1)
		{
		case 4://��һά��1,2,3,4�����
			for(int m=1; m<5; m++)
			{
				remainedTemp1.clear();
				removedTemp1.clear();
				remainedTemp1.push_back(m_commonFuncs.boundValue(m,4));
				remainedTemp1.push_back(m_commonFuncs.boundValue(m+1,4));
				removedTemp1.push_back(m_commonFuncs.boundValue(m+2,4));
				removedTemp1.push_back(m_commonFuncs.boundValue(m+3,4));
				remainedDimensionValues1.push_back(remainedTemp1);
				removedDimensionValues1.push_back(removedTemp1);
			}
			break;
		case 2://��һά��1,2��2,3��3,4��4,1�����
			remainedTemp1.clear();
			removedTemp1.clear();
			remainedTemp1.push_back(dimensionValues1.at(i).at(0));
			removedTemp1.push_back(dimensionValues1.at(i).at(1));
			remainedDimensionValues1.push_back(remainedTemp1);
			removedDimensionValues1.push_back(removedTemp1);

			remainedTemp1.clear();
			removedTemp1.clear();
			remainedTemp1.push_back(dimensionValues1.at(i).at(1));
			removedTemp1.push_back(dimensionValues1.at(i).at(0));
			remainedDimensionValues1.push_back(remainedTemp1);
			removedDimensionValues1.push_back(removedTemp1);
			break;
		case 1://��һάֻ��һ��ֵ�����ܽ�������
			break;
		}

		//���ÿһ�ֿ��ܳ��ֵ�����������µ�����Ϊ��������ȡԼ����ȥ����һ��������ĵ㣩
		for(int j=0; j<removedDimensionValues1.size(); j++)
		{
			//��ȡȥ������һ��ĵ㣬�����һ��ȥ���ж�������
			vector<CCell*> removedCells1;
			vector<vector<int>> tempDimensionValue1;

			//��ȡ��ȥ����һ��ĵ�
			m_commonFuncs.getTempDimensionValue(dimensionValues1,i,removedDimensionValues1.at(j),tempDimensionValue1);
			m_commonFuncs.ListAllCells(0,tempDimensionValue1,removedCells1);
			m_commonFuncs.kickOutCells(removedCells1,allOTPs);//ȥ����һ���в���ѡ��OTP
			/*
			*
			*ȥ��OTP֮������жϣ�
			*1)���removedCells�л��е㣬��˵����shrink��һ����Чshrink�������ж����ź͵ڶ����LIF��
			*2)���removedCellsΪ�ռ�������Ҫ����������ĵ��Լ��ڶ����LIF
			*/

			//��ʣ���
			if(removedCells1.size() > 0)
			{
				//m_commonFuncs.clearCellList(removedCells1);
				continue;
			}
			m_commonFuncs.clearCellList(removedCells1);
			//��ȡʣ�µ�һ���dimensionֵ
			vector<vector<int>> tempRemainedDimensionValue1;
			m_commonFuncs.getTempDimensionValue(dimensionValues1,i,remainedDimensionValues1.at(j),tempRemainedDimensionValue1);

			//�������ȥ�����г����п��ܵ�����
			vector<vector<CCell*>> expandCells1;
			for(int k=0; k<dimensionNum1; k++)//�����ж�k==i�����
			{	
				vector<CCell*> expandCell1;
				vector<CCell*> expandUTPCell;
				vector<vector<int>> tempExpandDimensionValue1;
				vector<int> expandDimensionValue1;
				int biggestValue1 = biggestDimensionValueAllowed1.at(k);//��ά��������ֵ��2��4��
				int kthRemainedDimensionValueNum1 = tempRemainedDimensionValue1.at(k).size();
				switch(kthRemainedDimensionValueNum1)
				{
				case 4://�����������һά������������
					break;
				case 2://����������������ʣ�µ�����������ԭ��Ϊ1,2�����䲿�־�Ϊ3,4
					if(k!=i&&biggestValue1==4)//1.k==iʱexpand֮���shrink֮ǰ��һ����;2.�����ά�����ȡֵֻ��2��Ҳ���ܽ���expand
					{
						int num1 = tempRemainedDimensionValue1.at(k).at(0);
						int num2 = tempRemainedDimensionValue1.at(k).at(1);
						if((num1==1&&num2==4) || (num1==4&&num2==1))
						{
							expandDimensionValue1.push_back(2);
							expandDimensionValue1.push_back(3);
						}
						else
						{
							int minNum = num1<num2? num1:num2;
							expandDimensionValue1.push_back(m_commonFuncs.boundValue(minNum+2,biggestValue1));
							expandDimensionValue1.push_back(m_commonFuncs.boundValue(minNum+3,biggestValue1));
						}
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,expandDimensionValue1,tempExpandDimensionValue1);
						m_commonFuncs.ListAllCells(0,tempExpandDimensionValue1,expandCell1);
						m_commonFuncs.keepCells(expandUTPCell ,expandCell1,allUTPs);
						m_commonFuncs.kickOutCells(expandCell1,allOTPs,allUTPs);//expand�ĵ㲻����TP
						if(expandCell1.size()>0){
							m_commonFuncs.clearCellList(expandCell1);
						}
						else{
							if (expandUTPCell.size()>0)
							{
								vector<vector<int>> dimensionValues2;
								vector<int> biggestDimensionValueAllowed2;
								//��ȡ��άֵ
								m_commonFuncs.getDimensionValue(dimensionValuesStr2,vOdd,dimensionValues2,biggestDimensionValueAllowed2);

								int dimensionNum2 = dimensionValues2.size();
								vector<vector<CCell*>> shrinkCells2;
								//��ÿһά������С����ȡ��С���ֵĵ�
								for(int kk=0; kk<dimensionNum2; kk++)
								{
									vector<CCell*> shrinkCell2;
									vector<vector<int>> tempShrinkDimensionValue2;
									vector<int> shrinkDimensionValue2;
									int ithDimensionNum2 = dimensionValues2.at(kk).size();
									switch(ithDimensionNum2)
									{
									case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
										for(int m=1; m<5; m++)
										{
											tempShrinkDimensionValue2.clear();
											shrinkCell2.clear();
											shrinkDimensionValue2.clear();
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
											m_commonFuncs.getTempDimensionValue(dimensionValues2,kk,shrinkDimensionValue2,tempShrinkDimensionValue2);
											m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue2,shrinkCell2);
											m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);
											m_commonFuncs.kickOutCommonCells(shrinkCell2,expandUTPCell);//Shrink�ĵ㲻���ǵ�һ��expand�ĵ�
											//m_commonFuncs.clearCellList(expandUTPCell);
											if(shrinkCell2.size()>0)
											{
												matrixCells.push_back(shrinkCell2);
											}
											else
											{
												m_commonFuncs.clearCellList(shrinkCell2);
											}
										}
										break;
									case 2: //�����������Ϊ1
										for(int m=0; m<2; m++)
										{
											tempShrinkDimensionValue2.clear();
											shrinkCell2.clear();
											shrinkDimensionValue2.clear();
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(kk).at(m),4));
											m_commonFuncs.getTempDimensionValue(dimensionValues2,kk,shrinkDimensionValue2,tempShrinkDimensionValue2);
											m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue2,shrinkCell2);
											m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);//Shrink�ĵ㲻����OTP
											m_commonFuncs.kickOutCommonCells(shrinkCell2,expandUTPCell);//Shrink�ĵ㲻���ǵ�һ��expand�ĵ�
											//m_commonFuncs.clearCellList(expandUTPCell);
											if(shrinkCell2.size()>0)
											{
												matrixCells.push_back(shrinkCell2);
											}
											else
											{
												m_commonFuncs.clearCellList(shrinkCell2);
											}
										}
										break;
									case 1://������������ٽ�������
										break;
									}
								}
							}
						}
					}
					m_commonFuncs.clearCellList(expandUTPCell);
					break;
				case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2��
					//���ǵÿ���k==iʱ��expand֮����ܺ�shrink֮ǰ��һ��
					int num = tempRemainedDimensionValue1.at(k).at(0);
					bool noPlus1 = false;
					bool noMinus1 = false;
					if(k==i)
					{

						if(biggestValue1==2)//���ʱ���ܽ���expand
						{
							noPlus1=true;
							noMinus1=true;
						}
						else
						{
							int anotherNum = removedDimensionValues1.at(j).at(0);
							if(m_commonFuncs.boundValue(num+1,biggestValue1)==anotherNum)
								noPlus1 = true;
							if(m_commonFuncs.boundValue(num-1,biggestValue1)==anotherNum)
								noMinus1 = true;
						}
					}
					else 
					{
						if(biggestValue1==2)//���ʱ��ֻ�ܵ�����expand
						{
							if(num==0)//�����ά���ֵ��2����ônumֻ������0��1
								noMinus1=true;
							else
								noPlus1=true;
						}
					}
					if(!noPlus1)
					{
						expandDimensionValue1.clear();
						tempExpandDimensionValue1.clear();
						expandCell1.clear();
						expandDimensionValue1.push_back(m_commonFuncs.boundValue(num+1,biggestValue1));
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,expandDimensionValue1,tempExpandDimensionValue1);
						m_commonFuncs.ListAllCells(0,tempExpandDimensionValue1,expandCell1);
						m_commonFuncs.keepCells(expandUTPCell ,expandCell1,allUTPs);
						m_commonFuncs.kickOutCells(expandCell1,allOTPs,allUTPs);//expand�ĵ㲻����TP
						if(expandCell1.size()>0){
							m_commonFuncs.clearCellList(expandCell1);
						}
						else{
							if (expandUTPCell.size()>0)
							{
								vector<vector<int>> dimensionValues2;
								vector<int> biggestDimensionValueAllowed2;
								//��ȡ��άֵ
								m_commonFuncs.getDimensionValue(dimensionValuesStr2,false,dimensionValues2,biggestDimensionValueAllowed2);

								int dimensionNum2 = dimensionValues2.size();
								vector<vector<CCell*>> shrinkCells2;
								//��ÿһά������С����ȡ��С���ֵĵ�
								for(int kk=0; kk<dimensionNum2; kk++)
								{
									vector<CCell*> shrinkCell2;
									vector<vector<int>> tempShrinkDimensionValue2;
									vector<int> shrinkDimensionValue2;
									int ithDimensionNum2 = dimensionValues2.at(kk).size();
									switch(ithDimensionNum2)
									{
									case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
										for(int m=1; m<5; m++)
										{
											tempShrinkDimensionValue2.clear();
											shrinkCell2.clear();
											shrinkDimensionValue2.clear();
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
											m_commonFuncs.getTempDimensionValue(dimensionValues2,kk,shrinkDimensionValue2,tempShrinkDimensionValue2);
											m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue2,shrinkCell2);
											m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);
											m_commonFuncs.kickOutCommonCells(shrinkCell2,expandUTPCell);//Shrink�ĵ㲻���ǵ�һ��expand�ĵ�
											//m_commonFuncs.clearCellList(expandUTPCell);
											if(shrinkCell2.size()>0)
											{
												matrixCells.push_back(shrinkCell2);
											}
											else
											{
												m_commonFuncs.clearCellList(shrinkCell2);
											}
										}
										break;
									case 2: //�����������Ϊ1
										for(int m=0; m<2; m++)
										{
											tempShrinkDimensionValue2.clear();
											shrinkCell2.clear();
											shrinkDimensionValue2.clear();
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(kk).at(m),4));
											m_commonFuncs.getTempDimensionValue(dimensionValues2,kk,shrinkDimensionValue2,tempShrinkDimensionValue2);
											m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue2,shrinkCell2);
											m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);//Shrink�ĵ㲻����OTP
											m_commonFuncs.kickOutCommonCells(shrinkCell2,expandUTPCell);//Shrink�ĵ㲻���ǵ�һ��expand�ĵ�
											//m_commonFuncs.clearCellList(expandUTPCell);
											if(shrinkCell2.size()>0)
											{
												matrixCells.push_back(shrinkCell2);
											}
											else
											{
												m_commonFuncs.clearCellList(shrinkCell2);
											}
										}
										break;
									case 1://������������ٽ�������
										break;
									}
								}
							}
						}
					}
					if(!noMinus1)
					{
						expandDimensionValue1.clear();
						tempExpandDimensionValue1.clear();
						expandCell1.clear();
						expandDimensionValue1.push_back(m_commonFuncs.boundValue(num-1,biggestValue1));
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,expandDimensionValue1,tempExpandDimensionValue1);
						m_commonFuncs.ListAllCells(0,tempExpandDimensionValue1,expandCell1);
						m_commonFuncs.keepCells(expandUTPCell ,expandCell1,allUTPs);
						m_commonFuncs.kickOutCells(expandCell1,allOTPs,allUTPs);//expand�ĵ㲻����TP
						if(expandCell1.size()>0){
							m_commonFuncs.clearCellList(expandCell1);
						}
						else{
							if (expandUTPCell.size()>0)
							{
								vector<vector<int>> dimensionValues2;
								vector<int> biggestDimensionValueAllowed2;
								//��ȡ��άֵ
								m_commonFuncs.getDimensionValue(dimensionValuesStr2,false,dimensionValues2,biggestDimensionValueAllowed2);

								int dimensionNum2 = dimensionValues2.size();
								vector<vector<CCell*>> shrinkCells2;
								//��ÿһά������С����ȡ��С���ֵĵ�
								for(int kk=0; kk<dimensionNum2; kk++)
								{
									vector<CCell*> shrinkCell2;
									vector<vector<int>> tempShrinkDimensionValue2;
									vector<int> shrinkDimensionValue2;
									int ithDimensionNum2 = dimensionValues2.at(kk).size();
									switch(ithDimensionNum2)
									{
									case 4: //�����������Ϊ2����1,2��2,3��3,4��4,1�������
										for(int m=1; m<5; m++)
										{
											tempShrinkDimensionValue2.clear();
											shrinkCell2.clear();
											shrinkDimensionValue2.clear();
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m,4));
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(m+1,4));
											m_commonFuncs.getTempDimensionValue(dimensionValues2,kk,shrinkDimensionValue2,tempShrinkDimensionValue2);
											m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue2,shrinkCell2);
											m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);
											m_commonFuncs.kickOutCommonCells(shrinkCell2,expandUTPCell);//Shrink�ĵ㲻���ǵ�һ��expand�ĵ�
											//m_commonFuncs.clearCellList(expandUTPCell);
											if(shrinkCell2.size()>0)
											{
												matrixCells.push_back(shrinkCell2);
											}
											else
											{
												m_commonFuncs.clearCellList(shrinkCell2);
											}
										}
										break;
									case 2: //�����������Ϊ1
										for(int m=0; m<2; m++)
										{
											tempShrinkDimensionValue2.clear();
											shrinkCell2.clear();
											shrinkDimensionValue2.clear();
											shrinkDimensionValue2.push_back(m_commonFuncs.boundValue(dimensionValues2.at(kk).at(m),4));
											m_commonFuncs.getTempDimensionValue(dimensionValues2,kk,shrinkDimensionValue2,tempShrinkDimensionValue2);
											m_commonFuncs.ListAllCells(0,tempShrinkDimensionValue2,shrinkCell2);
											m_commonFuncs.kickOutCells(shrinkCell2,allOTPs);//Shrink�ĵ㲻����OTP
											m_commonFuncs.kickOutCommonCells(shrinkCell2,expandUTPCell);//Shrink�ĵ㲻���ǵ�һ��expand�ĵ�
											//m_commonFuncs.clearCellList(expandUTPCell);
											if(shrinkCell2.size()>0)
											{
												matrixCells.push_back(shrinkCell2);
											}
											else
											{
												m_commonFuncs.clearCellList(shrinkCell2);
											}
										}
										break;
									case 1://������������ٽ�������
										break;
									}
								}
							}
						}
					}
					m_commonFuncs.clearCellList(expandUTPCell);
					break;
				}
			}
		}
	}

	return "";
}
