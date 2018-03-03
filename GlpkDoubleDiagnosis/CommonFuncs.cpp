#include "StdAfx.h"
#include "CommonFuncs.h"


ResourceGuard* CCommonFuncs:: rgd = nullptr ;
CCommonFuncs::CCommonFuncs(void)
{
	if(CCommonFuncs::rgd == nullptr)
	{
		CCommonFuncs::rgd = new ResourceGuard();
	}
}

CCommonFuncs::~CCommonFuncs(void)
{
}

//ʵ��nȡ1
void CCommonFuncs::chooseOne(vector<int> src, vector<int> &choosed)
{
	for(int i=0; i<src.size(); i++)
	{
		choosed.push_back(src.at(i));
	}
}

//�ݹ��оٳ���src��ȡx��Ԫ�ص����п��ܣ�ͨ������chooseOne����ʵ�֣���C5ȡ2�൱��C5ȡ1��C4ȡ1
void CCommonFuncs::chooseX(vector<int> src, int x, vector<vector<int>> &allChoosed)
{
	if(x>src.size())
		return;
	if(x==1)
	{
		vector<int> choosed;
		chooseOne(src,choosed);
		for(int i=0; i<choosed.size(); i++)
		{
			vector<int> v;
			v.push_back(choosed.at(i));
			allChoosed.push_back(v);
		}
	}
	else
	{
		vector<vector<int>> allChoosedTemp;
		chooseX(src,x-1,allChoosedTemp);

		vector<vector<int>> allRemainedTemp;
		getRemained(src,allChoosedTemp,allRemainedTemp);

		for(int i=0; i<allRemainedTemp.size(); i++)
		{
			vector<int> choosed;
			chooseOne(allRemainedTemp.at(i),choosed);

			for(int j=0; j<choosed.size(); j++)
			{
				vector<int> tempV;

				//�ȸ�����push
				allChoosed.push_back(allChoosedTemp.at(i));
				allChoosed.at(allChoosed.end()-1-allChoosed.begin()).push_back(choosed.at(j));

				/*for(int k=0; k<allChoosedTemp.at(i).size(); k++)
				{
				tempV.push_back(allChoosedTemp.at(i).at(k);
				}
				tempV.push_back(choosed.at(j);

				allChoosed.push_back(tempV);*/
			}
		}

		//ȥ�أ���Ϊ��˳���޹صģ�����123��213��231�ȶ���һ����
		for(int i=0; i<allChoosed.size(); i++)
		{
			int j=i+1;
			while(j<allChoosed.size())
			{	
				int k=0;
				int temp = allChoosed.at(j).size();//��j�е�Ԫ�ظ���
				while(k<temp)
				{
					if(find(allChoosed.at(i).begin(),allChoosed.at(i).end(),allChoosed.at(j).at(k))==allChoosed.at(i).end())
						break;
					k++;
				}
				if(k==temp)//��j���е�Ԫ���ڵ�i���ж��У����ظ���
					allChoosed.erase(allChoosed.begin()+j);
				else
					j++;
			}
		}
	}
}

//��ȡsrc�б�ȡʣ�µ�Ԫ��
void CCommonFuncs::getRemained(vector<int> src, vector<vector<int>> allChoosed, vector<vector<int>> &allRemained)
{
	for(int i=0; i<allChoosed.size(); i++)
	{
		vector<int> remained;
		for(int j=0; j<src.size(); j++)
		{
			if(find(allChoosed.at(i).begin(),allChoosed.at(i).end(),src.at(j))!=allChoosed.at(i).end())
				continue;
			remained.push_back(src.at(j));
		}
		allRemained.push_back(remained);
	}
}

//����choosed�е��У�ѡȡmatrix�е��Ӿ��󣬾����б�Ŵ�0��ʼ
void CCommonFuncs::getSubMatrix(vector<vector<int>> matrix, vector<int> choosed, vector<vector<int>> &subMatrix)
{
	if(choosed.size()==0)
		return;
	//�ȼ���choosed�еĵ�һ��
	int index = choosed.at(0);
	for(int i=0; i<matrix.size(); i++)
	{
		vector<int> row;
		row.push_back(matrix.at(i).at(index));
		subMatrix.push_back(row);
	}

	for(int i=1; i<choosed.size(); i++)
	{
		index = choosed.at(i);
		for(int j=0; j<matrix.size(); j++)
		{
			subMatrix.at(j).push_back(matrix.at(j).at(index));
		}
	}
}

//����choosed�е��У�ѡȡmatrix�е��Ӿ��󣬾����б�Ŵ�0��ʼ
void CCommonFuncs::getSubMatrix(vector<int*> matrix, vector<int> choosed, vector<vector<int>> &subMatrix)
{
	if(choosed.size()==0)
		return;
	//�ȼ���choosed�еĵ�һ��
	int index = choosed.at(0);
	for(int i=0; i<matrix.size(); i++)
	{
		vector<int> row;
		row.push_back(matrix.at(i)[index]);
		subMatrix.push_back(row);
	}

	for(int i=1; i<choosed.size(); i++)
	{
		index = choosed.at(i);
		for(int j=0; j<matrix.size(); j++)
		{
			subMatrix.at(j).push_back(matrix.at(j)[index]);
		}
	}

}

//��rowת��Ϊstring
string CCommonFuncs::RowTostring(int*row, int num)
{
	string rowStr="";
	for(int i=0; i<num; i++)
	{
		rowStr.append("%d ", row[i]);
	}
	return rowStr;
}

//���ַ����ָ�
vector<string> CCommonFuncs::split(string str, string pattern)
{
	std::string::size_type pos;
	std::vector<std::string> result;
	str += pattern;//��չ�ַ����Է������  
	int size = str.size();

	for (int i = 0; i<size; i++)
	{
		pos = str.find(pattern, i);
		if (pos<size)
		{
			std::string s = str.substr(i, pos - i);
			result.push_back(s);
			i = pos + pattern.size() - 1;
		}
	}
	return result;
}

//�ݹ���DimensionValue�б�ʾ������cell
void CCommonFuncs::ListAllCells(int k, vector<vector<int>> DimensionValue, vector<CCell*>& cellList)
{
	int dimension = DimensionValue.size();
	CCell* pCell;
	if(k==dimension-1)
	{
		for(int i=0; i<DimensionValue.at(k).size(); i++)
		{
			pCell = new CCell(dimension);
			CCommonFuncs::rgd->addResource(pCell);
			pCell->value[k] = DimensionValue.at(k).at(i);
			cellList.push_back(pCell);
		}
	}
	else
	{
		ListAllCells(k+1,DimensionValue,cellList);
		if(DimensionValue.at(k).size()>0)
		{
			for(int i=0; i<cellList.size(); i++)
			{
				cellList.at(i)->value[k] = DimensionValue.at(k).at(0);
			}
		}

		int cellNum = cellList.size();
		//��������cell
		for(int j=0; j<cellNum; j++)
		{
			for(int i=1; i<DimensionValue.at(k).size(); i++)
			{
				pCell = new CCell(dimension);
				CCommonFuncs::rgd->addResource(pCell);
				pCell->Copy(cellList.at(j));
				pCell->value[k] = DimensionValue.at(k).at(i);
				cellList.push_back(pCell);
			}
		}
	}
}

//��ȡÿһά�Ŀ���ȡֵ���Լ�tModel��ÿһά��������ֵ��һ�������ÿһά��������ֵΪ4�������������Ϊ����vOdd=true�����һά��������ֵΪ2��
void CCommonFuncs::getDimensionValue(string dimensionValuesStr, bool vOdd, vector<vector<int>>& dimensionValues, vector<int>& biggestDimensionValueAllowed)
{
	if(dimensionValuesStr.size()==0)
		return;
	//printf("enter getDimensionValue\n");
	vector<string> dimensions = split(dimensionValuesStr,"\n");
	int dimensionNum = dimensions.size();

	for(int i=0; i<dimensionNum; i++)
	{
		vector<int> dimension;
		vector<string> subStrs = split(dimensions[i]," ");
		for(int j=0; j<subStrs.size(); j++)
		{
			int temp = atoi(subStrs.at(j).c_str());
			if(temp>4 || temp<1)
			{
				//MessageBox("ÿά�Ŀ���ȡֵ��1��4�����������д��������䣡");
				if(temp != 0)dimensionValues.clear();
				return;
			}
			//printf("%d ", temp);
			dimension.push_back(temp);
		}
		//printf("\n");
		dimensionValues.push_back(dimension);

		if(vOdd && i==dimensionNum-1)
			biggestDimensionValueAllowed.push_back(2);
		else
			biggestDimensionValueAllowed.push_back(4);
	}
	//printf("out getDimensionValue\n");
}

//��value��ֵ�޶���tModel��ȡֵ��Χ�ڣ���1,2,3,4�е�һ��(1��4������һ���)
int CCommonFuncs::boundValue(int value, int biggestValue)
{
	while(value<0)
		value+=100;

	//biggestValueֻ����2��4
	int result = value%biggestValue;
	if(result==0)
		return biggestValue;
	else
		return result;
}

//��dimensionValue�еĵ�iά��ΪithTempDimensionValue���õ�����ֵ���浽 tempDimensionValue��
void CCommonFuncs::getTempDimensionValue(vector<vector<int>> dimensionValue, int i, vector<int> ithTempDimensionValue, 
	vector<vector<int>>& tempDimensionValue)
{
	int dimensionNum = dimensionValue.size();
	for(int m=0; m<dimensionNum; m++)
	{
		if(m!=i)
			tempDimensionValue.push_back(dimensionValue.at(m));
		else
			tempDimensionValue.push_back(ithTempDimensionValue);
	}
}

//��ȡmatrixCells�г��ֵ�����tModelIndex
void CCommonFuncs::getAlltModelIndexAppearedInMatrix(vector<vector<CCell*>> matrixCells, vector<int>& tModelIndexAppeared)
{
	for(int i=0; i<matrixCells.size(); i++)
	{
		for(int j=0; j<matrixCells.at(i).size(); j++)
		{
			int tModelIndex = matrixCells.at(i).at(j)->getTModelIndex();
			vector<int>::iterator iter = find(tModelIndexAppeared.begin(),tModelIndexAppeared.end(),tModelIndex);
			if(iter==tModelIndexAppeared.end())
				tModelIndexAppeared.push_back(tModelIndex);
		}
	}
}

//�������е��У���tModelIndexAppeared��,��matrixCellsת��Ϊϵ������
void CCommonFuncs::getCoefficientMatrix(vector<int> tModelIndexAppeared, vector<vector<CCell*>> matrixCells, vector<int*>& coefficientMatrix)
{
	for(int i=0; i<matrixCells.size(); i++)
	{
		int *coefficients = new int[tModelIndexAppeared.size()];
		CCommonFuncs::rgd->addResource(coefficients);
		for(int j=0; j<tModelIndexAppeared.size(); j++)//��ʼ��
			coefficients[j]=0;
		for(int j=0; j<matrixCells.at(i).size(); j++)//����Ӧ��ϵ��ֵ��Ϊ1
		{
			int tModelIndex = matrixCells.at(i).at(j)->getTModelIndex();
			vector<int>::iterator iter = find(tModelIndexAppeared.begin(),tModelIndexAppeared.end(),tModelIndex);
			int index = iter - tModelIndexAppeared.begin();
			coefficients[index]=1;
		}
		coefficientMatrix.push_back(coefficients);
	}
}

//�������е��У���tModelIndexAppeared��,��matrixCellsת��Ϊϵ������
void CCommonFuncs::getCoefficientMatrix(vector<int> tModelIndexAppeared, vector<vector<CCell*>> matrixCells, vector<vector<int>>& coefficientMatrix)
{
	for(int i=0; i<matrixCells.size(); i++)
	{
		vector<int> coefficients;
		for(int j=0; j<tModelIndexAppeared.size(); j++)//��ʼ��
			coefficients.push_back(0);
		for(int j=0; j<matrixCells[i].size(); j++)//����Ӧ��ϵ��ֵ��Ϊ1
		{
			int tModelIndex = matrixCells[i][j]->getTModelIndex();
			vector<int>::iterator iter = find(tModelIndexAppeared.begin(), tModelIndexAppeared.end(), tModelIndex);
			int index = iter - tModelIndexAppeared.begin();
			coefficients[index] = 1;
		}
		coefficientMatrix.push_back(coefficients);
	}
}

//��һάֵ��01��ʾ��ת��Ϊ1,2,3,4�ı�ʾ������11��3��ʾ
void CCommonFuncs::transDimensionLineFrom01totModelIndex(string ZeroOneLine, string& tModelIndexLine)
{
	vector<string> subStrs = split(ZeroOneLine," ");
	tModelIndexLine = "";
	int num = subStrs.size();
	for(int i=0; i<num; i++)
	{
		string tModelIndex;
		if(subStrs.at(i)=="00")
		{
			tModelIndex = "1";
		}
		else if(subStrs.at(i)=="01")
		{
			tModelIndex = "2";
		}
		else if(subStrs.at(i)=="11")
		{
			tModelIndex = "3";
		}
		else if(subStrs.at(i)=="10")
		{
			tModelIndex = "4";
		}
		else if(subStrs.at(i)=="0")
		{
			tModelIndex = "1";
		}
		else if(subStrs.at(i)=="1")
		{
			tModelIndex = "2";
		}

		tModelIndexLine.append(tModelIndex);
		if(i!=num-1)
			tModelIndexLine.append(" ");
	}
}
void CCommonFuncs::transDimensionLineFrom01totModelIndex(string ZeroOneLine, vector<int>& dimensionTModelIndex)
{
	vector<string> subStrs = split(ZeroOneLine," ");
	dimensionTModelIndex.clear();
	int num = subStrs.size();
	for(int i=0; i<num; i++)
	{
		int tModelIndex=0;
		if(subStrs.at(i)=="00")
		{
			tModelIndex = 1;
		}
		else if(subStrs.at(i)=="01")
		{
			tModelIndex = 2;
		}
		else if(subStrs.at(i)=="11")
		{
			tModelIndex = 3;
		}
		else if(subStrs.at(i)=="10")
		{
			tModelIndex = 4;
		}
		else if(subStrs.at(i)=="0")
		{
			tModelIndex = 1;
		}
		else if(subStrs.at(i)=="1")
		{
			tModelIndex = 2;
		}

		dimensionTModelIndex.push_back(tModelIndex);
	}
}


//�����01��ʾת��Ϊ1,2,3,4��ʾ�����110101��ʾ��322
int CCommonFuncs::getTModelIndexOf01Str(string ZeroOneStr, bool vOdd)
{
	int tModelIndex=0;

	//����ά����
	int strLength = ZeroOneStr.size();
	int dimensionNum = strLength/2;
	string tempStr;//���ں����ļ���
	if(vOdd)
		tempStr = ZeroOneStr.substr(0,strLength-1);
	else
		tempStr = ZeroOneStr;

	//��ά����ת���������һά��ʼ
	string subStr;
	int tempStrLength = tempStr.size();
	for(int i=0; i<dimensionNum; i++)
	{
		subStr = tempStr.substr(tempStrLength-2-i*2,2);
		if(subStr == "00")
			tModelIndex += 1*pow(10.0,i);
		else if(subStr == "01")
			tModelIndex += 2*pow(10.0,i);
		else if(subStr == "11")
			tModelIndex += 3*pow(10.0,i);
		else if(subStr == "10")
			tModelIndex += 4*pow(10.0,i);
	}

	//���Ϊ�����������һά
	if(vOdd)
	{
		tModelIndex *= 10;
		subStr = ZeroOneStr.substr(strLength-1,1);
		if(subStr == "0")
			tModelIndex+=1;
		else
			tModelIndex+=2;
	}

	return tModelIndex;
}

//����tModelIndexs��cellList�еĵ�ȥ��
void CCommonFuncs::kickOutCells(vector<CCell*>& cellList, vector<int> tModelIndexs)
{
	for(int i=0; i<cellList.size(); i++)
	{
		int tModelIndex = cellList.at(i)->getTModelIndex();
		if(find(tModelIndexs.begin(),tModelIndexs.end(),tModelIndex)!=tModelIndexs.end())//�����tModelIndexs����ȥ��
		{
			cellList.erase(cellList.begin()+i);
			i--;
		}
	}
}
//����tModelIndexs��cellList2�еĵ㱣�浽celllist1��
void CCommonFuncs::keepCells(vector<CCell*>& cellList1 ,vector<CCell*> cellList2,vector<int> tModelIndexs)
{
	for(int i=0; i<cellList2.size(); i++)
	{
		int tModelIndex = cellList2.at(i)->getTModelIndex();
		if(find(tModelIndexs.begin(),tModelIndexs.end(),tModelIndex)!=tModelIndexs.end())//�����tModelIndexs����ȥ��
		{
			cellList1.push_back(cellList2.at(i));
		}
	}
}

//����tModelIndexs��cellList1�еĵ㱣�浽cellList2��
void CCommonFuncs::getShrinkOtps(vector<CCell*>& cellList1, vector<CCell*>& cellList2, vector<int> tModelIndexs)
{
	for(int i=0; i<cellList1.size(); i++)
	{
		int tModelIndex = cellList1.at(i)->getTModelIndex();
		if(find(tModelIndexs.begin(),tModelIndexs.end(),tModelIndex)!=tModelIndexs.end())//�����tModelIndexs�������cellList2
		{
			cellList2.push_back(cellList1.at(i));
		}
	}

}

//����tModelIndexs1��tModelIndexs2��cellList�еĵ�ȥ��
void CCommonFuncs::kickOutCells(vector<CCell*>& cellList, vector<int> tModelIndexs1, vector<int> tModelIndexs2)
{
	for(int i=0; i<cellList.size(); i++)
	{
		int tModelIndex = cellList.at(i)->getTModelIndex();
		if(find(tModelIndexs1.begin(),tModelIndexs1.end(),tModelIndex)!=tModelIndexs1.end()||
			find(tModelIndexs2.begin(),tModelIndexs2.end(),tModelIndex)!=tModelIndexs2.end())//�����tModelIndexs1��tModelIndexs2����ȥ��
		{
			cellList.erase(cellList.begin()+i);
			i--;
		}
	}
}
//ȡcelllist1��cellList2�����ĵ㣬���浽cellList3��
void CCommonFuncs::getOverlapCells(vector<CCell*>& cellList1, vector<CCell*>& cellList2,vector<CCell*>& cellList3)
{
	for (int i = 0; i < cellList1.size(); i++)
	{
		for (int j = 0; j  < cellList2.size(); j ++)
		{
			if(cellList1.at(i)->compare(cellList2.at(j)))
			{
				cellList3.push_back(cellList1.at(i));
				break;
			}	
		}
	}
}
//��celllist1ɾ����cellList2�ظ��ĵ�
void CCommonFuncs::kickOutCommonCells(vector<CCell*>& cellList1, vector<CCell*> cellList2)
{
	for (int i = 0; i < cellList1.size(); i++)
	{
		for (int j = 0; j  < cellList2.size(); j ++)
		{
			if(cellList1.at(i)->compare(cellList2.at(j)))
			{
				cellList1.erase(cellList1.begin()+i);
				i--;
				break;
			}
			
		}

	}
}

//��1,2,3,4��ʾ�Ĵ�ת��Ϊ01��
string CCommonFuncs::transTModelIndexStrTo01Str(string tModelIndexStr, bool vOdd)
{
	string result = "";
	int len = tModelIndexStr.size();
	if(vOdd)
	{
		for(int i=0; i<len-1; i++)
		{
			if(tModelIndexStr.at(i)=='1')
				result.append("00");
			else if(tModelIndexStr.at(i)=='2')
				result.append("01");
			else if(tModelIndexStr.at(i)=='3')
				result.append("11");
			else if(tModelIndexStr.at(i)=='4')
				result.append("10");
		}

		//��������Ϊ����ʱ�����һά�����
		int i=len-1;
		if(tModelIndexStr.at(i)=='1')
			result.append("0");
		else if(tModelIndexStr.at(i)=='2')
			result.append("1");
	}
	else
	{
		for(int i=0; i<len; i++)
		{
			if(tModelIndexStr.at(i)=='1')
				result.append("00");
			else if(tModelIndexStr.at(i)=='2')
				result.append("01");
			else if(tModelIndexStr.at(i)=='3')
				result.append("11");
			else if(tModelIndexStr.at(i)=='4')
				result.append("10");
		}
	}
	return result;
}

void CCommonFuncs::transStrMatrix2VectorMatrix(string strMatrix, vector<vector<int>>& vectorMatrix)
{
	vector<string> lines = split(strMatrix,"\n");
	int lineNum = lines.size();
	for(int i=0; i<lineNum; i++)
	{
		vector<int> lineVector;
		vectorMatrix.push_back(lineVector);
	}

	for(int i=0; i<lineNum; i++)
	{
		vector<string> cols = split(lines.at(i)," ");
		for(int j=0; j<cols.size(); j++)
		{
			vectorMatrix.at(i).push_back(atoi(cols.at(j).c_str()));
		}
	}
}

//��vector��ʾ�ľ���ת��Ϊstring��ʽ
void CCommonFuncs::transVectorMatrix2StrMatrix(vector<vector<int>> vectorMatrix, string &strMatrix)
{
	strMatrix="";
	int lineNum = vectorMatrix.size();
	if(lineNum==0)
		return;
	int colNum = vectorMatrix.at(0).size();
	for(int i=0; i<lineNum; i++)
	{
		for(int j=0; j<colNum; j++)
		{
			string tempStr;
			tempStr.append("%d",vectorMatrix.at(i).at(j));
			strMatrix.append(tempStr);
			if(j==colNum-1)
				strMatrix.append("\r\n");
			else
				strMatrix.append(" ");
		}
	}
}

//�ж�cell�Ƿ���cellList��
bool CCommonFuncs::isCellInCellList(CCell* pCell, vector<CCell*> cellList)
{
	int num = cellList.size();
	for(int i=0; i<num; i++)
	{
		if(pCell->compare(cellList.at(i)))
			return true;
	}
	return false;
}

//���celllist
bool CCommonFuncs::clearCellList(vector<CCell*> &cellList)
{
	for (vector<CCell*>::iterator it = cellList.begin(); it != cellList.end(); it ++) 
	{
		if (NULL != *it&&(*it)->dimension>0) 
		{
			delete *it; 
			*it = NULL;
		}
	}
	cellList.clear();
	vector<CCell*>().swap(cellList);
	return true;
}