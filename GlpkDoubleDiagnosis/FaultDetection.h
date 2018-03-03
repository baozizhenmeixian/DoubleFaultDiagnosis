/*
������װ�˼��ȱ�ݵĽӿ�
*/
#pragma once
#include <vector>
using namespace std;
class FaultDetection
{
public:
	FaultDetection(void);
	~FaultDetection(void);
	//����:testset������������ԭ���ʽexprOri��������ʽexpr����Ч��testcase��Ӧ����validTCsol����Ч��testcase��Ӧ����invalidTCsol
	//����:true,������ʽû�з���ȱ�ݣ�false��������ʽ����ȱ��
	bool checkTestSet(vector<string> testset,string exprOri,string expr,vector<int>& validTCsol,vector<int>& invalidTCsol,bool vOdd);
	//���ݲ��������������ʽ�Ľ�����¾���
	//����:CoefficientMatrixԼ������TModelIndex�����е���������Ч�Ĳ�������validTCsol����Ч�Ĳ�������invalidTCsol
	//���أ�CoefficientMatrixԼ������ʣ�µ�����
	int updateMatrix
		(
		vector<vector<int>>& CoefficientMatrix,vector<int>& TModelIndex, 
		vector<int> validTCsol,vector<int> invalidTCsol,vector<int>& constrain_index,vector<bool>& mutantsremained);
	//����������ʽ�Ƿ�ȼ�
	bool checkEqualment(string str1, string str2);
	//ȥ��������ȫ����0������
	int removeZeroRow(vector<vector<int>> old_CoefficientMatrix,vector<vector<int>>& new_CoefficientMatrix);

	//ʹ���󽻼��ķ�ʽ������
	string checkExprwithSet(string exprOri, string expr, vector<string> testset, vector<vector<int>> case_mutant, vector<vector<int>>finalCoefficientMatrix, vector<int>finalTModelIndex, vector<string>& mutants_DP, bool vOdd);

	bool checkTestCase(int testcase_Index, string exprOri, string expr, bool vOdd);
private:
	//ԭ���ʽ
	string exprOri;
	//������ʽ
	string expr;

};
void removeAllOneCols(vector<vector<int>>& coefficientMatrix,vector<int> row_index);

