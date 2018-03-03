#include "StdAfx.h"
#include "DoubleLIRFHandler.h"


//����LRF�㷨
string CDoubleLIRFHandler::run2(
	string dimensionValuesStr, bool vOdd,
	vector<vector<CCell*>>& matrixCells,
	vector<int> allOTPs, vector<int> allUTPs
	)
{
	vector<int>literalNum;
	return runWithParameter2(0,literalNum, dimensionValuesStr,vOdd, matrixCells, allOTPs,  allUTPs);
}
//����LIF+LRF�㷨
string CDoubleLIRFHandler::runWithParameter2(
	int param, vector<int>&literalNum, 
	string dimensionValuesStr, bool vOdd, 
	vector<vector<CCell*>>& matrixCells, 
	vector<int> allOTPs, vector<int> allUTPs
	)
{
	int i=0;
	int	j=0;
	int k=0;
	int l=0;
	int n=0;
	int flag=0;

	vector<vector<int>> dimensionValues;
	vector<int> biggestDimensionValueAllowed;
	int allShrinksCount = 0;
	int validShrinksCount = 0;

	m_commonFuncs.getDimensionValue(dimensionValuesStr,vOdd,dimensionValues,biggestDimensionValueAllowed);

	int dimensionNum = dimensionValues.size();

	//���ÿһά����Ϊԭ����1/2����ѡһά��Ϊԭ����1/2
	for(i=0; i<dimensionNum; i++)
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

		//���ÿһ�ֿ��ܳ��ֵ�����������µ��ٴ�����
		for(j=0; j<removedDimensionValues.size(); j++)
		{
			//��ȡȥ������һ��ĵ㣬�����һ��ȥ���ж�������
			vector<CCell*> removedCells;
			vector<vector<int>> tempDimensionValue;
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
				matrixCells.push_back(matrixRow);*/
				continue;
			}

			//��ȡʣ�µ�һ���dimensionֵ
			vector<vector<int>> tempRemainedDimensionValue;
			m_commonFuncs.getTempDimensionValue(dimensionValues,i,remainedDimensionValues.at(j),tempRemainedDimensionValue);

			for(k=0; k<dimensionNum; k++)
			{
				vector<vector<int>> remainedDimensionValues1;
				vector<vector<int>> removedDimensionValues1;
				vector<int> remainedTemp1;
				vector<int> removedTemp1;

				int kthRemainedDimensionValueNum = tempRemainedDimensionValue.at(k).size();
				//�о���һά���ܳ��ֵ������������¼����
				switch(kthRemainedDimensionValueNum)
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
					remainedTemp1.push_back(tempRemainedDimensionValue.at(k).at(0));
					removedTemp1.push_back(tempRemainedDimensionValue.at(k).at(1));
					remainedDimensionValues1.push_back(remainedTemp1);
					removedDimensionValues1.push_back(removedTemp1);

					remainedTemp1.clear();
					removedTemp1.clear();
					remainedTemp1.push_back(tempRemainedDimensionValue.at(k).at(1));
					removedTemp1.push_back(tempRemainedDimensionValue.at(k).at(0));
					remainedDimensionValues1.push_back(remainedTemp1);
					removedDimensionValues1.push_back(removedTemp1);
					break;
				case 1://��һάֻ��һ��ֵ�����ܽ�������
					break;
				}

				for(l=0; l<removedDimensionValues1.size(); l++)
				{
					//��ȡȥ������һ��ĵ㣬�����һ��ȥ���ж�������
					vector<CCell*> removedCells1;
					vector<vector<int>> tempDimensionValue1;

					//��ȡ��ȥ����һ��ĵ�
					m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue,k,removedDimensionValues1.at(l),tempDimensionValue1);
					m_commonFuncs.ListAllCells(0,tempDimensionValue1,removedCells1);
					m_commonFuncs.kickOutCells(removedCells1,allOTPs);//ȥ����һ���в���ѡ��OTP

					/*
					*
					*ȥ��OTP֮������жϣ�
					*1)���removedCells�л��е㣬��˵����shrink��һ����Чshrink������������ĵ㡣
					*2)���removedCellsΪ�ռ�������Ҫ����������ĵ�
					*/
					//��ʣ���
					if(removedCells1.size() > 0)
					{
						//һ��
						vector<CCell*> matrixRow;
						for(int m=0; m<removedCells1.size(); m++)
						{
							matrixRow.push_back(removedCells1.at(m));
						}
						matrixCells.push_back(matrixRow);
						continue;
					}

					//��ȡʣ�µ�һ���dimensionֵ
					vector<vector<int>> tempRemainedDimensionValue1;
					m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue,k,remainedDimensionValues1.at(l),tempRemainedDimensionValue1);
					//�������ȥ�����г����п��ܵ�����
					vector<vector<CCell*>> expandCells;
					int literalN=0;
					for(n=0; n<dimensionNum; n++)
					{
						if(i==2&&j==2&&k==2&&l==1&&n==2){
							flag=1;

						}
						vector<CCell*> expandCell;
						vector<vector<int>> tempExpandDimensionValue;
						vector<int> expandDimensionValue;
						int biggestValue = biggestDimensionValueAllowed.at(n);//��ά��������ֵ��2��4��
						int kthRemainedDimensionValueNum1 = tempRemainedDimensionValue1.at(n).size();
						switch(kthRemainedDimensionValueNum1)
						{
						case 4://�����������һά������������
							break;
						case 2://����������������ʣ�µ�����������ԭ��Ϊ1,2�����䲿�־�Ϊ3,4
							if((n!=k&&biggestValue==4)/*||(n!=i&&biggestValue==4)*/)//1.n==kʱexpand֮���shrink֮ǰ��һ����;2.�����ά�����ȡֵֻ��2��Ҳ���ܽ���expand
							{
								int num1 = tempRemainedDimensionValue1.at(n).at(0);
								int num2 = tempRemainedDimensionValue1.at(n).at(1);
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
								m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,n,expandDimensionValue,tempExpandDimensionValue);
								m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell);
								m_commonFuncs.kickOutCells(expandCell,allOTPs,allUTPs);//expand�ĵ㲻����TP
								if(expandCell.size()>0){
									expandCells.push_back(expandCell);
									literalNum.push_back(literalN);
								}
								literalN++;//
							}
							break;
						case 1://��������£����������䷽ʽ������ԭ��Ϊ1�����䲿�ֿ�����4��2��
							//���ǵÿ���n==iʱ��expand֮����ܺ�shrink֮ǰ��һ��
							int num = tempRemainedDimensionValue1.at(n).at(0);
							bool noPlus1 = false;
							bool noMinus1 = false;
							if(n==k)
							{

								if(biggestValue==2)//���ʱ���ܽ���expand
								{
									noPlus1=true;
									noMinus1=true;
								}
								else
								{
									int anotherNum = removedDimensionValues1.at(l).at(0);
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
								m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,n,expandDimensionValue,tempExpandDimensionValue);
								m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell);
								m_commonFuncs.kickOutCells(expandCell,allOTPs,allUTPs);//expand�ĵ㲻����TP
								if(expandCell.size()>0){
									expandCells.push_back(expandCell);
									literalNum.push_back(literalN);
								}
								literalN++;//
							}
							if(!noMinus1)
							{
								expandDimensionValue.clear();
								tempExpandDimensionValue.clear();
								expandCell.clear();
								expandDimensionValue.push_back(m_commonFuncs.boundValue(num-1,biggestValue));
								m_commonFuncs.getTempDimensionValue(tempRemainedDimensionValue1,n,expandDimensionValue,tempExpandDimensionValue);
								m_commonFuncs.ListAllCells(0,tempExpandDimensionValue,expandCell);
								m_commonFuncs.kickOutCells(expandCell,allOTPs,allUTPs);//expand�ĵ㲻����TP
								if(expandCell.size()>0){
									expandCells.push_back(expandCell);
									literalNum.push_back(literalN);
								}
								literalN++;//
							}
							break;
						}
					}
					//�Ѿ���ȡ���е�expand cells����ȥ����һ��Ͷ�Ӧ�������һ��ϲ���ɾ����е�һ��
					for(int m=0; m<expandCells.size(); m++)
					{
						vector<CCell*> matrixRow;
						for(int a=0; a<removedCells.size(); a++)
							matrixRow.push_back(removedCells.at(a));
						for(int b=0; b<expandCells.at(m).size(); b++)
							matrixRow.push_back(expandCells.at(m).at(b));
						if (param!=1)
							matrixCells.push_back(matrixRow);
					}
				}
			}
		}
	}
	return "";
}

