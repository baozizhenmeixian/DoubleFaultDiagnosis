#include "stdafx.h"
#include "Utility.h"
#include "FaultDetection.h"
#include "CommonFuncs.h"
#include "BDD.h"
#include "GLPK_Solve.h"
#include "MatrixPartition.h"
#include<algorithm>
#include "public.h"

FaultDetection::FaultDetection(void)
{
}


FaultDetection::~FaultDetection(void)
{
}

string FaultDetection::checkExprwithSet(
	string exprOri/*ԭ���ʽ*/,
	string expr/*����*/,
	vector<string> testset/*��������*/,
	vector<vector<int>> case_mutant/*ע�⣬�������Լ������������ת��ϵ��λ�ñ�ţ�������һ�оͶ�Ӧһ���񸲸�*/,
	vector<vector<int>>finalCoefficientMatrix/*�����Ӧ��Լ�����󷽳�*/,
	vector<int>finalTModelIndex/*���в���������Ӧ�ĸ񸲸�*/,
	vector<string>& mutants_DP/*�����������ҵ���Ӧ�ı�����ʹ�õ�����*/,
	bool vOdd)
{
	GLPK_Solve glpk_Solve;
	CCommonFuncs commonFuncs;

	//1.ȱ�����ɲ�������
	vector<vector<int>> CoefficientMatrix(finalCoefficientMatrix);
	vector<int> TModelIndex(finalTModelIndex);/*�������в���������Ӧ�ĸ񸲸�*/
	//2.��������ʽ
	vector<int> validTCsol;
	vector<int> invalidTCsol;
	//�ò���������������ʽ���������Ч����û�з���ȱ��
	if (checkTestSet(testset/*��������*/, exprOri/*ԭ���ʽ*/, expr/*����*/, 
		validTCsol/*��һ�½ڵ㼯�ϣ�Ϊ��Ч���Ե�*/, invalidTCsol/*һ�½ڵ㼯�ϣ�Ϊ��Ч���Ե�*/, vOdd))
	{//�����棬˵��ȫ����Ч���Ե�
		return "no fault";
	}
	//�ò���������������ʽ��������Ч���������������ȱ��
	else
	{//���ؼ٣�˵������Ч���Ե�
		vector<int> p_mutants;
		vector<int> u_mutants;
		sort(invalidTCsol.begin(),invalidTCsol.end());
		vector<vector<int>> CoefficientMatrix_mid;
		vector<int> invalid_index;
		//vector<int> valid_index;
		//printf("enter FaultDetection While\n");
		while(true)
		{
			invalid_index.clear();
			/*printf("TModelIndex��\n");
			for (int j = 0; j < TModelIndex.size(); j++)
			{
				printf("%d ", TModelIndex[j]);
			}
			printf("\nend TModelIndex\r\n");
			printf("invalidTCsol��\n");
			for (int j = 0; j < invalidTCsol.size(); j++)
			{
				printf("%d ", invalidTCsol[j]);
			}
			printf("\nend invalidTCsol\r\n");
			printf("case_mutant��\n");
			for (int j = 0; j < case_mutant.size(); j++)
			{
				for (int l = 0; l < case_mutant[j].size(); l++)
				{
					printf("%d ", case_mutant[j][l]);
				}
				printf("\n");
			}
			printf("\nend case_mutant\r\n");*/
			//�󽻼�
			for(int i=0;i<invalidTCsol.size();)
			{
				for(int j=0;j<TModelIndex.size()&&i<invalidTCsol.size();j++)
				{
					if(invalidTCsol[i]==TModelIndex[j])
					{//�ڼ������ҵ������ϵ�ֵ�����Ĳ��������ڵ�
						invalid_index.push_back(j);//�����в����������ҵ���Ӧ�񸲸�֮���¼����
						if(p_mutants.empty())
						{
							p_mutants=case_mutant[j];//��ȡӵ�иø񸲸ǵ�����ϵ������ı��
						}
						else
						{
							vector<int> t_mutants2;
							set_intersection(/*�󽻼�*/
								p_mutants.begin(),/*���м������*/
								p_mutants.end(),/*���м����յ�*/
								case_mutant[j].begin(),/*�����Ӧ��Լ������ϵ�����*/
								case_mutant[j].end(),/*�����Ӧ��Լ������ϵ���յ�*/
								inserter(
									t_mutants2, 
									t_mutants2.end()
								)/*��װ��std���뼯�Ϻ�������һ������Ϊ���Ͽռ䣬�ڶ�������Ϊ����λ��*/
							);//��װ��std�󼯺Ͻ�������
							p_mutants.swap(t_mutants2);/*�����󵽼����뼯�Ͻ���*/
							vector<int>().swap(t_mutants2);/*�����ʱ����*/
						}
						if(p_mutants.size()==1)
							break;//�������ֻ��һ����¼���˳�������Ϊʲô�أ���Ϊ�����ζ���ҵ���Ψһ�ı���Լ������
						i++;
					}
				}
				if(p_mutants.size()==1)
					break;
			}//��ν�󽻼������ǰ������д���񸲸ǵ�Լ������������Ϊ��ת��
			printf("p_mutants1��\n");//p_mutants���Ƿ��̵ı��
			for (int j = 0; j < p_mutants.size(); j++)
			{
				printf("%d ", p_mutants[j]);
			}
			printf("\nend p_mutants1\r\n");
			printf("�󽻼�����\n");
			//�󲢼�
			/*printf("validTCsol��\n");
			for (int j = 0; j < validTCsol.size(); j++)
			{
				printf("%d ", validTCsol[j]);
			}
			printf("\nend validTCsol\r\n");*/
			for(int i=0;i<validTCsol.size()&&p_mutants.size()>1;)
			{
				for(int j=0;j<TModelIndex.size()&&i<validTCsol.size();j++)
				{
					if(validTCsol[i]==TModelIndex[j])
					{
						//valid_index.push_back(j);
						vector<int> t_mutants;
						set_difference(//��
							p_mutants.begin(),
							p_mutants.end(),
							case_mutant[j].begin(),
							case_mutant[j].end(),
							inserter(
								t_mutants, 
								t_mutants.end()
							)
						);
						p_mutants.swap(t_mutants);
						vector<int>().swap(t_mutants); 
						if(p_mutants.size()==1)
							break;
						i++;
					}
				}
				if(p_mutants.size()==1)
					break;
			}//��ν������ʵ����˵������󽻼��õ��Ľ����Ψһ����Ҫ����Щ����"һ�½ڵ�"��Լ�������޳�
			printf("p_mutants2��\n");
			for (int j = 0; j < p_mutants.size(); j++)
			{
				printf("%d ", p_mutants[j]);
			}
			printf("\nend p_mutants2\r\n");
			printf("�󲢼�����\n");
			//�ɴ˿����ƶ�����������߼�
			//���Ȼ�ȡ����㣨Ҳ�����������ʽ��һ�µĵ㣩�ͷǴ���㣨Ҳ�����������ʽһ�µĵ㣩
			//������Լ������������в�һ�µķ��̹鲢������A��ʹ�õ�����ת���������A�а�����һ�µĵ㣬��Ҫ�޳�
			//A+(S-B)���ǿ�����Ҫ�жϵļ��ϵ㣬�´ε�����ʱ���ʹ�ô˼�����Ϊȫ��S������ͬ�Ĵ���
			if(p_mutants.size()==1)
				break;

			//���´������˵��������ǵò���Ψһ��Լ�����̣������²�����������
			//���ɾ���
			CoefficientMatrix_mid.clear(); 
			removeAllOneCols(CoefficientMatrix,p_mutants);
			for(int i=0;i<p_mutants.size();i++)
			{//p_mutants���ھͳ��˴��ų����̵��������
				bool zero=true;
				for(int j=0;j<CoefficientMatrix[p_mutants[i]].size();j++)
				{//CoefficientMatrix����Լ����������CoefficientMatrix[p_mutants[i]]���Ǵ���Լ������
					vector<int>::iterator it = find( invalid_index.begin(), invalid_index.end(), j); //����Լ�������и����ڵ��Ƿ����Ѿ�У����ĵ�
					if(it!=invalid_index.end())
					{//�����ȣ����ʾû�в鵽������ȣ����ʾ�鵽�ˣ�˵�����������Ѿ�У����ĵ㣬�޸�Լ�������иýڵ��ֵΪ0�����ڱ�ʾ�޳��˸ýڵ㣨�񸲸ǣ�
						CoefficientMatrix[p_mutants[i]][j]=0;
					}
					if(CoefficientMatrix[p_mutants[i]][j]!=0)
					{//û�в鵽���߲���У����ĵ�
						zero=false;
						/*bool same=true;
						for(int ii=0;ii<p_mutants.size();ii++)
						{
							if(CoefficientMatrix[ii][j]==0)
							{
								same=false;
							}
						}*///���ô���
					}

				}
				if(!zero)
				{
					CoefficientMatrix_mid.push_back(CoefficientMatrix[p_mutants[i]]);//���µõ���Լ�����̲��뵽һ����ʱ�ռ���
				}
			}
			printf("���ɾ������\n"); 
			if(CoefficientMatrix_mid.empty())
				break;
			//2.2.���ɲ�������
			vector<vector<int>> CoefficientMatrix_temp;
			vector<int> index_temp;
			removeAllZeroCols(CoefficientMatrix_mid, TModelIndex, CoefficientMatrix_temp, index_temp);
			CoefficientMatrix_mid.clear();
			removeZeroRow(CoefficientMatrix_temp, CoefficientMatrix_mid);
			if(CoefficientMatrix_mid.empty())
				break;
			if(CoefficientMatrix_mid.size()==p_mutants.size()&&index_temp.size()==1)
				break;
			else
			{
				printf("����������û�еõ�Ψһ���\n");
				//testset = testCaseObj.genTestCaseWithGLPK(CoefficientMatrix_mid, index_temp);
				//�����������ε�һ��
				vector<int> colSelected;
				glpk_Solve.run(colSelected, &CoefficientMatrix_mid);//���½������ʱ��Լ������
				//���õ��ĸ񸲸�ת���ɲ�������
				//printf("tModelIndexAppeared colSelected Print:\n");
				testset.clear();
				for (int j = 0; j < colSelected.size(); j++)
				{
					char tmp[40];
					sprintf(tmp, "%d", index_temp[colSelected[j]]);
					string inde = "";
					inde.append(tmp);
					string point01Str = commonFuncs.transTModelIndexStrTo01Str(inde, vOdd);
					testset.push_back(point01Str);
					//printf("%s\n", point01Str.c_str());
				}

				printf("\n");
				//�����������ε�һ�����
				//2.3.��������ʽ
				validTCsol.clear();
				invalidTCsol.clear();
				checkTestSet(testset,exprOri,expr,validTCsol,invalidTCsol,vOdd);//��������
			}
			printf("���ɲ�����������\n");
		}
		//����ʣ�µ�һ�ж�Ӧ�ı���
		if (p_mutants.size() == 0){//���һ����û�У���ֱ�ӷ���һ������������
			return mutants_DP[0];
		}
		int mutant_index=p_mutants[0];
		printf("%s %d mutant_index:%d\n", __FILE__, __LINE__, mutant_index);
		printf("�õ��Ľ��[%d]:%s\n", mutant_index, mutants_DP[mutant_index].c_str());
		return mutants_DP[mutant_index];
	}
}
//��������߼���ʲô�أ�
bool FaultDetection::checkTestSet(
	vector<string> testset/*���ɵĲ����������ϣ������������ȫ��ͨ�������ʾû������*/, 
	string exprOri/*ԭ���ʽ*/, 
	string expr/*����*/, 
	vector<int>& validTCsol, /*�����һ�µĽڵ㣬Ϊ��Ч���Ե�*/
	vector<int>& invalidTCsol/*���һ�µĽڵ㣬Ϊ��Ч���Ե�*/, 
	bool vOdd)
{
	Utility uti;
	CCommonFuncs CCommonFuncs;

	string oriexp = exprOri;
	string faultexp = expr;
	printf("checkTestSet testset��\n");
	for (int i = 0; i < testset.size();i++){
		printf("testset[%d]:%s\n",i,testset[i].c_str());
	}
	printf("checkTestSet testset end\n");

	printf("exprOri:%s\n", exprOri.c_str());
	printf("expr:%s\n", expr.c_str());
	for(int i=0;i<testset.size();i++)
	{
		string testcase = testset[i];
		printf("testcase:%s\n", testcase.c_str());
		printf("oriexp testcase evaluate:%d\n", uti.evaluate(oriexp, testcase));//���������ڴ˱��ʽ���Ƿ�Ϊ��
		printf("oriexp faultexp testcase evaluateMutant:%d\n", uti.evaluateMutant(oriexp, faultexp, testcase));//���������ڴ˱��ʽ���Ƿ�Ϊ��
		bool v=(uti.evaluate(oriexp,testcase))^(uti.evaluateMutant(oriexp,faultexp,testcase));//����ֵ��ͬ��Ϊ�٣���ͬ��Ϊ��,������������Դ���˷�ʱ��ѽ�������и��ŵ��㷨
		printf("v:%d\n", v);
		if(v){
			invalidTCsol.push_back(CCommonFuncs.getTModelIndexOf01Str(testset[i],vOdd));
		}
		else
		{
			validTCsol.push_back(CCommonFuncs.getTModelIndexOf01Str(testset[i],vOdd));
		}
	}
	if(invalidTCsol.size()==0)
		return true;
	else
		return false;
}
int FaultDetection::updateMatrix
	(
	vector<vector<int>>& CoefficientMatrix,vector<int>& TModelIndex, 
	vector<int> validTCsol,vector<int> invalidTCsol,vector<int>& constrain_index,vector<bool>& mutantsremained)
{
	vector<vector<int>> temp_CoefficientMatrix;
	vector<vector<int>> old_CoefficientMatrix(CoefficientMatrix);
	vector<int> temp_TModelIndex;
	vector<bool> temp_mutantsremained;
	vector<int> temp_constrain_index;
	vector<bool> allzero_cow;
	for(int i=0;i<CoefficientMatrix.size()&&CoefficientMatrix.size()>0;i++)
	{
		bool allzero = true;
		for(int j=0;j<CoefficientMatrix[i].size();j++)
		{
			if(CoefficientMatrix[i][j]==1)
			{
				//�ж�������ǲ���valid�ĵ�
				vector<int>::iterator valid = find( validTCsol.begin(), validTCsol.end(), TModelIndex[j]); 	
				//����Ч�ĵ㣬���ᷢ�����ȱ��
				if(valid!=validTCsol.end())
				{
					mutantsremained[constrain_index[i]]=false;//��ʾ���ᷢ�����ȱ��
					break;
				}
				//�ж�������ǲ�����Ч�ĵ�
				vector<int>::iterator invalid = find( invalidTCsol.begin(), invalidTCsol.end(), TModelIndex[j]); 
				//����Ч�ĵ㣬��ϵ����Ϊ0
				if (invalid!=invalidTCsol.end( ))
				{
					if(mutantsremained[constrain_index[i]])
						mutantsremained[constrain_index[i]]=true;//��ʾ���ܷ������ȱ��
					for(int k=0;k<CoefficientMatrix.size();k++)
					{
						if(old_CoefficientMatrix[k][j]==0)
						{
							mutantsremained[constrain_index[k]]=false;
						}
					}
					CoefficientMatrix[i][j]=0;
				}
			}
			else
			{
				//�ж�������ǲ�����Ч�ĵ�
				vector<int>::iterator invalid = find( invalidTCsol.begin(), invalidTCsol.end(), TModelIndex[j]); 
				//����Ч�ĵ㣬��ϵ����Ϊ0
				if (invalid!=invalidTCsol.end( ))
				{
					mutantsremained[constrain_index[i]]=false;//��ʾ�����ܷ������ȱ��
				}
			}
		}
		for(int j=0;j<CoefficientMatrix[i].size();j++)
		{
			if(CoefficientMatrix[i][j]==1)
			{
				allzero=false;
			}
		}
		allzero_cow.push_back(allzero);
	}
	for(int i=0;i<CoefficientMatrix.size();i++)
	{
		if(mutantsremained[constrain_index[i]]/*&&!allzero_cow[i]*/)
		{
			temp_CoefficientMatrix.push_back(CoefficientMatrix[i]);
			temp_constrain_index.push_back(constrain_index[i]);
		}
	}

	constrain_index.swap(temp_constrain_index);
	temp_TModelIndex = TModelIndex;
	CoefficientMatrix.clear();
	TModelIndex.clear();
	for(int i=0; i<temp_CoefficientMatrix.size(); i++)
	{
		vector<int> coefficient;
		CoefficientMatrix.push_back(coefficient);
	}

	for(int i=0; i<temp_TModelIndex.size(); i++)
	{
		//�鿴���е�ÿһ���ǲ��Ƕ�Ϊ0
		int j=0;
		int lineNum = temp_CoefficientMatrix.size();
		for(; j<lineNum; j++)
		{
			if(temp_CoefficientMatrix.at(j).at(i)==1)
				break;
		}
		if(j<lineNum)//���в�Ϊȫ0���ӵ�finalMatrix��
		{
			for(int m=0; m<lineNum; m++)
				CoefficientMatrix.at(m).push_back(temp_CoefficientMatrix.at(m).at(i));

			TModelIndex.push_back(temp_TModelIndex.at(i));
		}
	}
	int mutantcount = 0;
	for (int i = 0; i < mutantsremained.size(); i++)
	{
		if(mutantsremained[i])
		{
			mutantcount+=1;
		}
	}
	return mutantcount;
}
bool FaultDetection::checkEqualment(string str1, string str2)
{
	BDD string1(str1);
	BDD string2(str2);

	return string1.isEquivalent(string2);
}

int FaultDetection::removeZeroRow(vector<vector<int>> old_CoefficientMatrix,vector<vector<int>>& new_CoefficientMatrix)
{
	
	for(int i=0; i<old_CoefficientMatrix.size();i++)
	{
		bool allzero = true;
		for (int j= 0; j < old_CoefficientMatrix[i].size(); j++)
		{
			if(old_CoefficientMatrix[i][j]==1)
			{
				allzero=false;
				break;
			}
		}
		if(!allzero)
		{
			new_CoefficientMatrix.push_back(old_CoefficientMatrix[i]);
		}
	}
	return new_CoefficientMatrix.size();
}

void removeAllOneCols(vector<vector<int>>& coefficientMatrix,vector<int> row_index)
{
	vector<int> allOneCol;
	bool have_valid_col=false;
	for(int i=0; i<coefficientMatrix[0].size(); i++)
	{
		//�鿴���е�ÿһ���ǲ��Ƕ�Ϊ0
		int j=0;
		int lineNum = row_index.size();
		int oneCount=0;
		for(; j<lineNum; j++)
		{
			if(coefficientMatrix.at(row_index[j]).at(i)==1)
				oneCount++;
		}
		if(oneCount==lineNum)//����ȫΪ1
		{
			allOneCol.push_back(i);
		}
		else
		{
			if(oneCount>0)
			{
				have_valid_col=true;
			}
		}
	}

	int col=have_valid_col?0:1;
	for(;col<allOneCol.size();col++)
	{
		for(int m=0;m<row_index.size();m++)
		{
			coefficientMatrix[row_index[m]][allOneCol[col]]=0;
		}
	}
}

bool FaultDetection::checkTestCase(int testcase_Index, string exprOri, string expr, bool vOdd)
{
	Utility uti;

	string oriexp = exprOri;
	string faultexp = expr;
	string testcasetemp;
	testcasetemp.append("%d",testcase_Index); 

	string testcase = testcasetemp;

	bool result=(uti.evaluate(oriexp,testcase))^(uti.evaluateMutant(oriexp,faultexp,testcase));

	return result;
}