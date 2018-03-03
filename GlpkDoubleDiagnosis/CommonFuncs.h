#pragma once
/* 2012.3 By �ν�ɼ
* �������ڴ��һЩ���õķ���
*/
#include <vector>
#include <algorithm>
#include "Cell.h"
#include <set>
using namespace std;

class
	ResourceGuard
{
public :
	void addResource(CCell* ccell_pointer){
		resources_cells.insert(ccell_pointer);

	}
	void addResource(int* int_array){
		resources_int_arrays.insert(int_array);

	}
	~ResourceGuard()
	{
		for ( std::set<CCell*>::iterator it = resources_cells.begin(); it!= resources_cells.end();it++)
		{
			delete (*it);
		}

		for ( std::set<int*>::iterator it = resources_int_arrays.begin(); it!= resources_int_arrays.end();it++)
		{
			delete[] (*it);
		}
	
	}
private:
	std::set<CCell*> resources_cells;
	std::set<int*> resources_int_arrays;
   
};
class CCommonFuncs
{
public:
	CCommonFuncs(void);
	~CCommonFuncs(void);
	//ʵ��nȡ1
	void chooseOne(vector<int> src, vector<int> &choosed);
	//�ݹ��оٳ���src��ȡx��Ԫ�ص����п��ܣ�ͨ������chooseOne����ʵ�֣���C5ȡ2�൱��C5ȡ1��C4ȡ1
	void chooseX(vector<int> src, int x, vector<vector<int>> &allChoosed);
	//��ȡsrc�б�ȡʣ�µ�Ԫ��
	void getRemained(vector<int> src, vector<vector<int>> allChoosed, vector<vector<int>> &allRemained);
	//����choosed�е��У�ѡȡmatrix�е��Ӿ��󣬾����б�Ŵ�0��ʼ
	void getSubMatrix(vector<vector<int>> matrix, vector<int> choosed, vector<vector<int>> &subMatrix);
	//��rowת��Ϊstring
	string RowTostring(int*row, int num);
	//��split�ָ�str
	vector<string> split(string str, string split);

	//�ݹ���DimensionValue�б�ʾ������cell
	void ListAllCells(int k, vector<vector<int>> DimensionValue, vector<CCell*>& cellList);

	//��ȡÿһά�Ŀ���ȡֵ���Լ�tModel��ÿһά��������ֵ��һ�������ÿһά��������ֵΪ4�������������Ϊ����vOdd=true�����һά��������ֵΪ2��
	void getDimensionValue(string dimensionValuesStr, bool vOdd, vector<vector<int>>& dimensionValues, 
		vector<int>& biggestDimensionValueAllowed);

	//��value��ֵ�޶���tModel��ȡֵ��Χ�ڣ���1,2,3,4�е�һ��(1��4������һ���)
	int boundValue(int value, int biggestValue);

	//��dimensionValue�еĵ�iά��ΪithTempDimensionValue���õ�����ֵ���浽 tempDimensionValue��
	void getTempDimensionValue(vector<vector<int>> dimensionValue, int i, vector<int> ithTempDimensionValue, 
		vector<vector<int>>& tempDimensionValue);

	//��ȡmatrixCells�г��ֵ�����tModelIndex
	void getAlltModelIndexAppearedInMatrix(vector<vector<CCell*>> matrixCells, 
		vector<int>& tModelIndexAppeared);

	//�������е��У���tModelIndexAppeared��,��matrixCellsת��Ϊϵ������
	void getCoefficientMatrix(vector<int> tModelIndexAppeared, vector<vector<CCell*>> matrixCells, 
		vector<int*>& coefficientMatrix);

	//�������е��У���tModelIndexAppeared��,��matrixCellsת��Ϊϵ������
	void getCoefficientMatrix(vector<int> tModelIndexAppeared, vector<vector<CCell*>> matrixCells, 
		vector<vector<int>>& coefficientMatrix);

	//��һάֵ��01��ʾ��ת��Ϊ1,2,3,4�ı�ʾ������11��3��ʾ
	void transDimensionLineFrom01totModelIndex(string ZeroOneLine, string& tModelIndexLine);
	void transDimensionLineFrom01totModelIndex(string ZeroOneLine, vector<int>& tModelIndexLine);

	//�����01��ʾת��Ϊ1,2,3,4��ʾ�����110101��ʾ��322
	int getTModelIndexOf01Str(string ZeroOneStr, bool vOdd);

	//����choosed�е��У�ѡȡmatrix�е��Ӿ��󣬾����б�Ŵ�0��ʼ
	void getSubMatrix(vector<int*> matrix, vector<int> choosed, vector<vector<int>> &subMatrix);

	//����tModelIndexs��cellList�еĵ�ȥ��
	void kickOutCells(vector<CCell*>& cellList, vector<int> tModelIndexs);
	//����tModelIndexs��cellList2�еĵ㱣�浽celllist1��
	void keepCells(vector<CCell*>& cellList1 ,vector<CCell*> cellList2,vector<int> tModelIndexs);

	//����tModelIndexs��cellList1�еĵ㱣�浽cellList2��
	void getShrinkOtps(vector<CCell*>& cellList1, vector<CCell*>& cellList2, vector<int> tModelIndexs);
	//ȡcelllist1��cellList2�����ĵ㣬���浽cellList3��
	void getOverlapCells(vector<CCell*>& cellList1, vector<CCell*>& cellList2,vector<CCell*>& cellList3);
	//��celllist1ɾ����cellList2�ظ��ĵ�
	void kickOutCommonCells(vector<CCell*>& cellList1, vector<CCell*> cellList2);

	//����tModelIndexs1��tModelIndexs2��cellList�еĵ�ȥ��
	void kickOutCells(vector<CCell*>& cellList, vector<int> tModelIndexs1, vector<int> tModelIndexs2);

	//��1,2,3,4��ʾ�Ĵ�ת��Ϊ01��
	string transTModelIndexStrTo01Str(string tModelIndexStr, bool vOdd);

	//������ʾ�ľ���ת��Ϊvector��ʾ�ľ���
	void transStrMatrix2VectorMatrix(string strMatrix, vector<vector<int>>& vectorMatrix);

	//��vector��ʾ�ľ���ת��Ϊstring��ʽ
	void transVectorMatrix2StrMatrix(vector<vector<int>> vectorMatrix, string &strMatrix);

	//�ж�cell�Ƿ���cellList��
	bool isCellInCellList(CCell* pCell, vector<CCell*> cellList);
	//���celllist
	bool clearCellList(vector<CCell*> &cellList);
	
	static ResourceGuard* rgd ;
};
