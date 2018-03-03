#include "StdAfx.h"
#include "DoubleLORFHandler.h"


//����LRF�㷨
string CDoubleLORFHandler::run2(
	string dimensionValuesStr, bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<int> allOTPs, vector<int> allUTPs
	)
{
	vector<int>literalNum;
	return runWithParameter2(0,literalNum, dimensionValuesStr,vOdd, matrixCells, allOTPs,  allUTPs);
}
//����LOF+LRF�㷨
string CDoubleLORFHandler::runWithParameter2(
	int param, vector<int>&literalNum, 
	string dimensionValuesStr, bool vOdd, 
	vector<vector<CCell*>>& matrixCells, 
	vector<int> allOTPs, vector<int> allUTPs
	)
{
	vector<vector<int>> dimensionValues;  //ԭ����
	vector<int> biggestDimensionValueAllowed;
	m_commonFuncs.getDimensionValue(dimensionValuesStr,vOdd,dimensionValues,biggestDimensionValueAllowed);

	int dimensionNum = dimensionValues.size();
	//���ÿһά����Ϊԭ����1/2����ѡһά����Ϊԭ��������

	for(int i=0; i<dimensionNum; i++) // ѡȡһά�ȣ���ʼ����
	{
		//�ֱ��ȡʣ�µ�һ���ȥ����һ��
		vector<vector<int>> remainedDimensionValues;
		vector<vector<int>> removedDimensionValues;
		vector<int> remainedTemp;
		vector<int> removedTemp;

		int ithDimensionValueNum = dimensionValues.at(i).size();
		//�о���һά���ܳ��ֵ������������¼����
		switch(ithDimensionValueNum)
		{
		case 4://��һά��1,2,3,4�����
			for(int m=1; m<5; m++)
			{
				remainedTemp.clear();
				removedTemp.clear();
				remainedTemp.push_back(m_commonFuncs.boundValue(m,4));
				remainedTemp.push_back(m_commonFuncs.boundValue(m+1,4));
				removedTemp.push_back(m_commonFuncs.boundValue(m+2,4));
				removedTemp.push_back(m_commonFuncs.boundValue(m+3,4));
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
		for(int j=0; j<removedDimensionValues.size(); j++)
		{
			//��ȡȥ������һ��ĵ㣬�����һ��ȥ���ж�������
			vector<CCell*> removedCells;
			vector<vector<int>> tempDimensionValue;//��һ���������ȡ�������õľ���

			//��ȡ��ȥ����һ��ĵ�
			m_commonFuncs.getTempDimensionValue(dimensionValues,i,removedDimensionValues.at(j),tempDimensionValue);
			m_commonFuncs.ListAllCells(0,tempDimensionValue,removedCells);
			m_commonFuncs.kickOutCells(removedCells,allOTPs);//ȥ����һ���в���ѡ��OTP
			/*
			*
			*ȥ��OTP֮������жϣ�
			*1)���removedCells�л��е㣬��˵����shrink��һ����Чshrink������������ĵ㡣
			*2)���removedCellsΪ�ռ�������Ҫ����������ĵ�
			*/

			//��ʣ���
			if(removedCells.size() > 0)
			{//��LIORF�Ĳ��������м��
				//һ��
				/*vector<CCell*> matrixRow;
				for(int m=0; m<removedCells.size(); m++)
				{
					matrixRow.push_back(removedCells.at(m));
				}
				if(param!=2)
					matrixCells.push_back(matrixRow);*/
				continue;
			}

			//��ȡʣ�µ�һ���dimensionֵ
			vector<vector<int>> tempRemainedDimensionValue1;//��һ�����������µľ���
			m_commonFuncs.getTempDimensionValue(dimensionValues,i,remainedDimensionValues.at(j),tempRemainedDimensionValue1);

			//�������ȥ�����г����п��ܵ�����
			vector<vector<CCell*>> expandCells1;
			vector<vector<CCell*>> expandCells2;
			vector<vector<int>> remainedDimensionValues1;//�洢���ź�����
			int literalN=0;//
			for(int k=0; k<dimensionNum; k++)//�����ж�k==i�����
			{	
				remainedDimensionValues1.clear();
				vector<int> remainedTemp1;
				vector<CCell*> expandCell;
				vector<vector<int>> tempExpandDimensionValue;
				vector<int> expandDimensionValue;
				int biggestValue = biggestDimensionValueAllowed.at(k);//��ά��������ֵ��2��4��
				int kthRemainedDimensionValueNum = tempRemainedDimensionValue1.at(k).size();
				switch(kthRemainedDimensionValueNum)
				{
				case 4://�����������һά������������
					break;
				case 2://����������������ʣ�µ�����������ԭ��Ϊ1,2�����䲿�־�Ϊ3,4
					if(k!=i&&biggestValue==4)//1.k==iʱexpand֮���shrink֮ǰ��һ����;2.�����ά�����ȡֵֻ��2��Ҳ���ܽ���expand
					{
						int num1 = tempRemainedDimensionValue1.at(k).at(0);
						int num2 = tempRemainedDimensionValue1.at(k).at(1);
						if((num1==1&&num2==4) || (num1==4&&num2==1))
						{
							expandDimensionValue.push_back(2);
							expandDimensionValue.push_back(3);
						}
						else
						{
							int minNum = num1<num2? num1:num2;
							expandDimensionValue.push_back(m_commonFuncs.boundValue(minNum+2,biggestValue));
							expandDimensionValue.push_back(m_commonFuncs.boundValue(minNum+3,biggestValue));
						}
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,expandDimensionValue,tempExpandDimensionValue);
						m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell);
						m_commonFuncs.kickOutCells(expandCell,allOTPs,allUTPs);//expand�ĵ㲻����TP
						if(removedCells.size()+expandCell.size()>0){
							expandCells1.push_back(expandCell);
							literalNum.push_back(literalN);
						}
						else
						{
							remainedTemp1.clear();
							remainedTemp1.push_back(1);
							remainedTemp1.push_back(2);
							remainedTemp1.push_back(3);
							remainedTemp1.push_back(4);
							remainedDimensionValues1.push_back(remainedTemp1);
						}
						literalN++;//
					}
					break;
				case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2��
					//���ǵÿ���k==iʱ��expand֮����ܺ�shrink֮ǰ��һ��
					int num = tempRemainedDimensionValue1.at(k).at(0);
					bool noPlus1 = false;
					bool noMinus1 = false;
					if(k==i)
					{

						if(biggestValue==2)//���ʱ���ܽ���expand
						{
							noPlus1=true;
							noMinus1=true;
						}
						else
						{
							int anotherNum = removedDimensionValues.at(j).at(0);
							if(m_commonFuncs.boundValue(num+1,biggestValue)==anotherNum)
								noPlus1 = true;
							if(m_commonFuncs.boundValue(num-1,biggestValue)==anotherNum)
								noMinus1 = true;
						}
					}
					else 
					{
						if(biggestValue==2)//���ʱ��ֻ�ܵ�����expand
						{
							if(num==0)//�����ά���ֵ��2����ônumֻ������0��1
								noMinus1=true;
							else
								noPlus1=true;
						}
					}
					if(!noPlus1)
					{
						expandDimensionValue.clear();
						tempExpandDimensionValue.clear();
						expandCell.clear();
						expandDimensionValue.push_back(m_commonFuncs.boundValue(num+1,biggestValue));
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,expandDimensionValue,tempExpandDimensionValue);
						m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell);
						m_commonFuncs.kickOutCells(expandCell,allOTPs,allUTPs);//expand�ĵ㲻����TP
						if(removedCells.size()+expandCell.size()>0){
							expandCells1.push_back(expandCell);
							literalNum.push_back(literalN);
						}
						else
						{
							//�洢���µ�
							remainedTemp1.clear();
							remainedTemp1.push_back(num);
							remainedTemp1.push_back(m_commonFuncs.boundValue(num+1,biggestValue));
							remainedDimensionValues1.push_back(remainedTemp1);
							//
						}
						literalN++;//
					}
					if(!noMinus1)
					{

						expandDimensionValue.clear();
						tempExpandDimensionValue.clear();
						expandCell.clear();
						expandDimensionValue.push_back(m_commonFuncs.boundValue(num-1,biggestValue));
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,expandDimensionValue,tempExpandDimensionValue);
						m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell);
						m_commonFuncs.kickOutCells(expandCell,allOTPs,allUTPs);//expand�ĵ㲻����TP
						if(removedCells.size()+expandCell.size()>0){
							expandCells1.push_back(expandCell);
							literalNum.push_back(literalN);
						}
						else
						{
							//�洢���µ�
							remainedTemp1.clear();
							remainedTemp1.push_back(m_commonFuncs.boundValue(num-1,biggestValue));
							remainedTemp1.push_back(num);
							remainedDimensionValues1.push_back(remainedTemp1);
							//
						}
						literalN++;//
					}
					break;
				}// ��һ�����Ž���
				if(remainedDimensionValues1.size()>0)
				{
					//��ʼ�ڶ�������
					for(int l=0;l<remainedDimensionValues1.size();l++)
					{
						//�����ź��dimensionֵ
						vector<vector<int>> tempRemainedDimensionValue2;//��һ�����ź�ľ���
						m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,k,remainedDimensionValues1.at(l),tempRemainedDimensionValue2);



						for(int n=0; n<dimensionNum; n++)//�����ж�k==i�����
						{	
							vector<int> remainedTemp1;
							vector<CCell*> expandCell1;
							vector<vector<int>> tempExpandDimensionValue;
							vector<int> expandDimensionValue;
							int biggestValue = biggestDimensionValueAllowed.at(n);//��ά��������ֵ��2��4��
							int kthRemainedDimensionValueNum = tempRemainedDimensionValue2.at(n).size();
							switch(kthRemainedDimensionValueNum)
							{
							case 4://�����������һά������������
								break;
							case 2://����������������ʣ�µ�����������ԭ��Ϊ1,2�����䲿�־�Ϊ3,4
								if(n!=i&&biggestValue==4)//1.k==iʱexpand֮���shrink֮ǰ��һ����;2.�����ά�����ȡֵֻ��2��Ҳ���ܽ���expand
								{
									int num1 = tempRemainedDimensionValue2.at(n).at(0);
									int num2 = tempRemainedDimensionValue2.at(n).at(1);
									expandDimensionValue.clear();
									tempExpandDimensionValue.clear();
									expandCell1.clear();
									if((num1==1&&num2==4) || (num1==4&&num2==1))
									{
										expandDimensionValue.push_back(2);
										expandDimensionValue.push_back(3);
									}
									else
									{
										int minNum = num1<num2? num1:num2;
										expandDimensionValue.push_back(m_commonFuncs.boundValue(minNum+2,biggestValue));
										expandDimensionValue.push_back(m_commonFuncs.boundValue(minNum+3,biggestValue));
									}
									m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue2,n,expandDimensionValue,tempExpandDimensionValue);
									m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell1);
									m_commonFuncs.kickOutCells(expandCell1,allOTPs,allUTPs);//expand�ĵ㲻����TP
									if(removedCells.size()+expandCell1.size()>0){
										expandCells2.push_back(expandCell1);
										literalNum.push_back(literalN);
									}
									literalN++;//
								}
								break;
							case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2��
								//���ǵÿ���k==iʱ��expand֮����ܺ�shrink֮ǰ��һ��
								int num = tempRemainedDimensionValue2.at(n).at(0);
								bool noPlus1 = false;
								bool noMinus1 = false;
								if(n==i)
								{

									if(biggestValue==2)//���ʱ���ܽ���expand
									{
										noPlus1=true;
										noMinus1=true;
									}
									else
									{
										int anotherNum = removedDimensionValues.at(j).at(0);
										if(m_commonFuncs.boundValue(num+1,biggestValue)==anotherNum)
											noPlus1 = true;
										if(m_commonFuncs.boundValue(num-1,biggestValue)==anotherNum)
											noMinus1 = true;
									}
								}
								else 
								{
									if(biggestValue==2)//���ʱ��ֻ�ܵ�����expand
									{
										if(num==0)//�����ά���ֵ��2����ônumֻ������0��1
											noMinus1=true;
										else
											noPlus1=true;
									}
								}
								if(!noPlus1)
								{
									expandDimensionValue.clear();
									tempExpandDimensionValue.clear();
									expandCell1.clear();
									expandDimensionValue.push_back(m_commonFuncs.boundValue(num+1,biggestValue));
									m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue2,n,expandDimensionValue,tempExpandDimensionValue);
									m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell1);
									m_commonFuncs.kickOutCells(expandCell1,allOTPs,allUTPs);//expand�ĵ㲻����TP
									if(removedCells.size()+expandCell1.size()>0){
										expandCells2.push_back(expandCell1);
										literalNum.push_back(literalN);
									}
									literalN++;//
								}
								if(!noMinus1)
								{

									expandDimensionValue.clear();
									tempExpandDimensionValue.clear();
									expandCell1.clear();
									expandDimensionValue.push_back(m_commonFuncs.boundValue(num-1,biggestValue));
									m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue2,n,expandDimensionValue,tempExpandDimensionValue);
									m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell1);
									m_commonFuncs.kickOutCells(expandCell1,allOTPs,allUTPs);//expand�ĵ㲻����TP
									if(removedCells.size()+expandCell1.size()>0){
										expandCells2.push_back(expandCell1);
										literalNum.push_back(literalN);
									}
									literalN++;//
								}
								break;
							}
						}
					}
				}
			}

			//�Ѿ���ȡ���е�expand cells����ȥ����һ��Ͷ�Ӧ�������һ��ϲ���ɾ����е�һ��
			for(int m=0; m<expandCells1.size(); m++)
			{
				vector<CCell*> matrixRow;
				for(int a=0; a<removedCells.size(); a++)
					matrixRow.push_back(removedCells.at(a));
				for(int b=0; b<expandCells1.at(m).size(); b++)
					matrixRow.push_back(expandCells1.at(m).at(b));
				if (param!=1)
					matrixCells.push_back(matrixRow);
			}
			for(int m=0; m<expandCells2.size(); m++)
			{
				if (param!=1)
					matrixCells.push_back(expandCells2.at(m));
			}
		}
	}//��һ��ѡȡһά����������


	return "";
}

