// BoolCheckFramework.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include <stdio.h>
#include "Expression.h"
#include "Mutants.h"
#include "DnfMutant.h"
#include "CommonFuncs.h"
#include <time.h>
#include "GLPK_Solve.h"
#include "LIFGenerateCellsForMutants.h"
#include "LOFGenerateCellsForMutants.h"
#include "LRFGenerateCellsForMutants.h"
#include "TOFGenerateCellsForMutants.h"
#include "Utility.h"
#include "DnfMutant.h"
#include "FaultDetection.h"
#include "public.h"
#include "DoubleLIRFHandler.h"
#include "DoubleLORFHandler.h"
#include "DoubleLRRFHandler.h"
#include "LRFHandler.h"
#include "DTLIFhandler.h"
#include "DTLIRFHandler.h"
#include "DTLRFHandler.h"
#include "DTTOFLRFhandler.h"
#include "DTDFCommonFuncs.h"

#define DNF 1
#define CNF 0

#pragma warning(disable:4996)

//���ݱ��ʽ��ȡ����
//����exp��Ϊԭ���ʽ
//����faultType��Ϊȱ������
//����iModel��Ϊ��ʽ����
//����DNF���ʽ
string generateDNF(vector<string> terms){
	vector<string>::iterator iter;
	for (int i = 0; i < terms.size(); i++){
		terms[i] = generateIDNForICNF(terms[i], true);
	}
	string oriExp;
	sort(terms.begin(), terms.end());
	iter = unique(terms.begin(), terms.end());
	if (iter != terms.end()){
		terms.erase(iter, terms.end());
	}
	for (int i = 0; i < terms.size(); i++){
		if (i == 0)oriExp = terms[i];
		else oriExp += "+" + terms[i];
	}
	return oriExp;
}



bool checkTestSet(
	vector<string> testset/*���ɵĲ����������ϣ������������ȫ��ͨ�������ʾû������*/,
	string exprOri/*ԭ���ʽ*/,
	string expr/*����*/,
	vector<int>& unformityTCsol, /*���һ�µĽڵ㣬Ϊ��Ч���Ե�*/
	vector<int>& differTCsol, /*�����һ�µĽڵ㣬Ϊ��Ч���Ե�*/
	vector<string>& unformityStr, /*���һ�µĽڵ㣬Ϊ��Ч���Ե�*/
	vector<string>& differTStr, /*�����һ�µĽڵ㣬Ϊ��Ч���Ե�*/
	vector<string>& unformityType,
	vector<string>& differType,
	bool vOdd)
{
	Utility uti;
	CCommonFuncs CCommonFuncs;

	string oriexp = exprOri;
	string faultexp = expr;
	for (int i = 0; i<testset.size(); i++)
	{
		string testcase = testset[i];
		bool v1 = (uti.evaluate(oriexp, testcase));
		bool v2 = (uti.evaluateMutant(oriexp, faultexp, testcase));
		bool v = v1 ^ v2;//����ֵ��ͬ��Ϊ�٣���ͬ��Ϊ��,������������Դ���˷�ʱ��ѽ�������и��ŵ��㷨
		if (v){
			if (v1 && !v2) differType.push_back("s");
			else differType.push_back("e");
			differTStr.push_back(testset[i]);
			differTCsol.push_back(CCommonFuncs.getTModelIndexOf01Str(testset[i], vOdd));
		}
		else
		{
			if (v1) unformityType.push_back("t");
			else unformityType.push_back("f");
			unformityStr.push_back(testset[i]);
			unformityTCsol.push_back(CCommonFuncs.getTModelIndexOf01Str(testset[i], vOdd));
		}
	}
	if (differTCsol.size() == 0)
		return true;
	else
		return false;
}

//�޳�Լ��������ȫ�����
void removeAllZeroRows(vector<vector<int>> finalCoefficientMatrix, vector<vector<int>>& coefficientMatrix, vector<int>& deletedRows){
	for (int i = 0; i < finalCoefficientMatrix.size(); i++){
		int l = 0;
		for (int j = 0; j < finalCoefficientMatrix[i].size(); j++){
			if (finalCoefficientMatrix[i][j] == 1){
				l = 1;
			}
		}
		if (l != 0){
			coefficientMatrix.push_back(finalCoefficientMatrix[i]);
		}
		else deletedRows.push_back(i);
	}
}

bool getLIFxLIFCells2(CCommonFuncs commonFuncs, vector<string> terms, int dimensionNum, vector<vector<string>>termsDimension, vector<string> opover2, vector<int> allOPs, bool vOdd, vector<vector<CCell *>>&matrixCells, vector<string>& mutants, vector<string>& faultTypes, vector<string>& faultTerms){
	if (terms.size() <= 1)
	{
		printf("���ʽֻ����һ������ܷ�������Double Fault LIFxLIF��");
		return false;
	}
	//double fault lifxlif ÿ����Ҫ��������
	// 1.1Ԥ������ÿ��term��dimensionValuesStr������
	vector<string> dimensionValuesStrs; //�洢ÿ��term��dimension value string
	int termNum = terms.size();
	for (int i = 0; i<termNum; i++)//���term���н���
	{
		string dimensionValuesStr = "";
		for (int j = 0; j < dimensionNum; j++)
		{
			//����ά��01��ʾ����ת��Ϊ1,2,3,4��ʾ��
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			if (j != dimensionNum - 1)
				dimensionValuesStr.append("\n");
		}
		dimensionValuesStrs.push_back(dimensionValuesStr);
	}
	//1.2 Ԥ�����ȡ��������otp
	vector<int> all3OPs;
	for (int i = 0; i<opover2.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(opover2.at(i), vOdd);
		all3OPs.push_back(value);
	}
	sort(all3OPs.begin(), all3OPs.end());

	for (int i = 0; i<termNum - 1; i++)
	{
		for (int j = i + 1; j < termNum; j++)
		{
			DTLIFhandler DTLIFhandler;
			int ln = matrixCells.size();
			DTLIFhandler.run2(dimensionValuesStrs.at(i), dimensionValuesStrs.at(j),
				matrixCells, allOPs, all3OPs);
			for (int k = ln; k < matrixCells.size(); k++){
				mutants.push_back("");
				faultTypes.push_back("LIFxLIF");
				faultTerms.push_back(terms[i] + "+" + terms[j]);
			}
		}
	}
	return true;
}
bool getLIFxLRFCells2(CCommonFuncs commonFuncs, vector<string> terms, int dimensionNum, vector<vector<string>>termsDimension, vector<string> opover2, vector<int> allUPs, vector<int> allOPs, bool vOdd, vector<vector<CCell *>>&matrixCells, vector<string>& mutants, vector<string>& faultTypes, vector<string>& faultTerms){
	if (terms.size() <= 1)
	{
		printf("���ʽֻ����һ������ܷ�������Double Fault LIFxLRF��");
		return false;
	}
	//double fault lifxlrf ÿ����Ҫ��������
	// 1.1Ԥ������ÿ��term��dimensionValuesStr������
	int termNum = terms.size();
	vector<string> dimensionValuesStrs; //�洢ÿ��term��dimension value string
	for (int i = 0; i<termNum; i++)//���term���н���
	{
		string dimensionValuesStr = "";
		for (int j = 0; j < dimensionNum; j++)
		{
			//����ά��01��ʾ����ת��Ϊ1,2,3,4��ʾ��
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			if (j != dimensionNum - 1)
				dimensionValuesStr.append("\n");
		}
		dimensionValuesStrs.push_back(dimensionValuesStr);
	}
	//1.2 Ԥ�����ȡ��������otp
	vector<int> all3OPs;
	for (int i = 0; i<opover2.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(opover2.at(i), vOdd);
		all3OPs.push_back(value);
	}
	sort(all3OPs.begin(), all3OPs.end());

	//����LRF����
	/*for (int i = 0; i<termNum; i++)
	{
		CLRFHandler lrfHandler;
		lrfHandler.run(dimensionValuesStrs.at(i), vOdd, matrixCells, allOPs, allUPs);
	}*/
	//
	for (int i = 0; i < termNum; i++)
	{
		bool hasEqualLRF = false;
		for (int j = 0; j < termNum; j++)
		{
			if (i == j)
			{
				continue;
			}
			CDTLIRFHandler dtlirfHandler;
			int ln = matrixCells.size();
			dtlirfHandler.run2(dimensionValuesStrs.at(i), dimensionValuesStrs.at(j), vOdd, matrixCells, allOPs, allUPs, hasEqualLRF);
			for (int k = ln; k < matrixCells.size(); k++){
				mutants.push_back("");
				faultTypes.push_back("LIFxLRF");
				faultTerms.push_back(terms[i] + "+" + terms[j]);
			}
		}
	}
	return true;
}
bool getLRFxLRFCells2(CCommonFuncs commonFuncs, vector<string> terms, int dimensionNum, vector<vector<string>>termsDimension, vector<int> allUPs, vector<int> allOPs, bool vOdd, vector<vector<CCell *>>&matrixCells, vector<string>& mutants, vector<string>& faultTypes, vector<string>& faultTerms){
	int termNum = terms.size();
	if (termNum <= 1)
	{
		printf("���ʽֻ����һ������ܷ�������Double Fault LRFxLRF��");
		return false;
	}
	//double fault lrfxlrf ÿ����Ҫ��������
	// 1.1Ԥ������ÿ��term��dimensionValuesStr������
	vector<string> dimensionValuesStrs; //�洢ÿ��term��dimension value string
	for (int i = 0; i<termNum; i++)//���term���н���
	{
		string dimensionValuesStr = "";
		bool isSingleLitral = false;
		//����ռȫ��ά����Ŀ
		int fullDemNum = 0;
		//����ռ���ά����Ŀ
		int halfDemNum = 0;
		for (int j = 0; j < dimensionNum; j++)
		{
			//����ά��01��ʾ����ת��Ϊ1,2,3,4��ʾ��
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			int strLength = tempStr.size();
			if (vOdd != true || j != dimensionNum - 1) //����ռȫʱ���ַ�������ӦΪ7,ռ�볤��Ϊ3
			{
				if (strLength == 7)
					fullDemNum++;
				else if (strLength == 3)
					halfDemNum++;
			}
			else //����ռȫʱ���ַ�������ӦΪ3,ռ�볤��Ϊ1
			{
				if (strLength == 3)
					fullDemNum++;
				else if (strLength == 1)
					halfDemNum++;
			}
			if (j != dimensionNum - 1)
				dimensionValuesStr.append("\n");
		}
		if (fullDemNum == dimensionNum - 1 && halfDemNum == 1)
			isSingleLitral = true;
		dimensionValuesStrs.push_back(dimensionValuesStr);
		//����DLRF����
		CDLRFHandler dlrfHandler;
		int ln = matrixCells.size();
		dlrfHandler.run2(dimensionValuesStr, vOdd, matrixCells, allOPs, allUPs);
		for (int k = ln; k < matrixCells.size(); k++){
			mutants.push_back("");
			faultTypes.push_back("LRFxLRF");
			faultTerms.push_back(terms[i]);
		}
	}
	return true;
}
string getLiterInTerm(string term, int literNo){
	for (int i = 0, j = 0; i < term.size() && j <= literNo; i++){
		if (term[i] == '!'){
			continue;
		}
		else{
			if (j == literNo){
				if (i > 0 && term[i - 1] == '!'){
					return string("!") + term[i];
				}
				else{
					return string({ term[i] });
				}
			}
			else j++;
		}
	}
	return "";
}
bool checkExistLiterInTerm(string term, string Liter){
	const char * ms = term.c_str();
	char ldt;
	if (Liter[0] == '!')ldt = Liter[1];
	else ldt = Liter[0];
	for (int i = 0; i < term.size(); i++){
		if (ms[i] == ldt){
			return true;
		}
	}
	return false;
}
int getDimensionNo(string model, string Liter, int &literalNo){
	const char * ms = model.c_str();
	int dimensionNo = -1;
	char ldt;
	if (Liter[0] == '!')ldt = Liter[1];
	else ldt = Liter[0];
	for (int i = 0; i < model.size(); i++){
		if (ms[i] == ldt){
			literalNo = i;
			dimensionNo = i / 2;
			break;
		}
	}
	return dimensionNo;
}
int getTermNo(vector<string> terms, string Term){
	int dimensionNo = -1;
	for (int i = 0; i < terms.size(); i++){
		if (terms[i] == Term){
			dimensionNo = i;
			break;
		}
	}
	return dimensionNo;
}
vector<int> checkShrinkTerms(string oriExp, vector<string>terms, vector<string>testcase){
	vector<int> faultTermsNo;
	Utility uti;
	for (int i = 0; i < terms.size(); i++){
		for (int j = 0; j < testcase.size(); j++){
			bool v1 = uti.evaluate(oriExp, testcase[j]);
			if (!v1)continue;//���������㣬��ȥ
			bool v2 = uti.evaluateMutant(oriExp, terms[i], testcase[j]);
			bool v = v1 ^ v2;//���ֵͬΪ�٣�ֵ��ͬΪ��
			if (v == false){//���Ϊ�٣�˵����������ȱ��
				faultTermsNo.push_back(i);
				break;
			}
		}
	}
	return faultTermsNo;
}
vector<vector<vector<string>>> changeDimensionToVector(vector<vector<string>> termsDimension){
	vector<vector<vector<string>>> a;
	CCommonFuncs commonFunc;
	for (int x = 0; x < termsDimension.size(); x++){
		vector<vector<string>> b;
		for (int y = 0; y < termsDimension[x].size(); y++){
			vector<string> c;
			c = commonFunc.split(termsDimension[x][y], " ");
			b.push_back(c);
			c.clear();
		}
		a.push_back(b);
		b.clear();
	}
	return a;
}
vector<vector<string>> getLOFChangeCells(string expression, string model, int termNo, string Liter1){
	//�ֽ���ʽ
	PRE_PROCESS_DATA preData;

	//ת��ΪIDNF
	Utility uti;
	vector<string> terms = uti.extractDNFstyle(expression);
	vector<string>::iterator iter;
	string newExp = "";
	//�õ�IDNF
	newExp = expression;// generateIDNForICNF(expression, true);
	//printf("newExp:%s\n", newExp.c_str());

	getPreProcessData(newExp, preData, DNF);
	CCommonFuncs commonFuncs;
	int literalNum = model.size();
	int dimensionNum = (literalNum + 1) / 2;
	bool vOdd = (literalNum % 2 == 0) ? false : true;
	//printf("%d,%d,%d\n", __LINE__, vOdd, literalNum);

	//��ȡ�仯�ʵ�ά�ȱ��
	int literalNo1;
	int dimensionNo1 = getDimensionNo(model, Liter1, literalNo1);

	//��ȡ���е�OTP
	vector<string> op = preData.op;
	vector<int> allOPs;
	for (int i = 0; i<op.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(op.at(i), vOdd);
		allOPs.push_back(value);
	}
	sort(allOPs.begin(), allOPs.end());

	//��ȡ���е�UTP
	vector<int> allUPs;
	vector<string> up = preData.up;
	for (int i = 0; i<up.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(up.at(i), vOdd);
		allUPs.push_back(value);
	}
	sort(allUPs.begin(), allUPs.end());

	//��ȡ�Ƿ��е��ʵ��ж�
	bool isSingleLitral = false;

	/*printf("preData.termsDimension:\n");
	for (int x = 0; x < preData.termsDimension.size(); x++){
		for (int y = 0; y < preData.termsDimension[x].size(); y++){
			printf("%s\n", preData.termsDimension[x][y].c_str());
		}
		printf("\n\n");
	}
	printf("\n");*/

	vector<vector<vector<string>>>termsDimensions = changeDimensionToVector(preData.termsDimension);
	vector<vector<string>> termsDimension = termsDimensions[termNo];
	string mutantTerm = terms[termNo];
	if (!checkExistLiterInTerm(mutantTerm, Liter1)){
		printf("cannot occur lof!\n");
		vector<vector<string>> a;
		return a;
	}
	if (mutantTerm.size() == 1 || mutantTerm.size() == 2 && mutantTerm[0] == '!'){
		printf("just can occur tof");
		vector<vector<string>> a;
		return a;
	}
	string mutantChangeTerm;
	//һ��term��Ӧ��ά��ֵ���м���ά�ȣ�dimensionNum�ͻ��ж��
	vector<string> testcase;
	CCommonFuncs commonFunc;
	mutantChangeTerm = mutantTerm;
	string lm = Liter1;
	int lz = mutantChangeTerm.find(lm);
	string::iterator itr = mutantChangeTerm.begin() + lz;
	mutantChangeTerm.erase(itr);
	//printf("LOF:mutantChangeTerm[%s]\n", mutantChangeTerm.c_str());

	vector<string> addValues;
	//�����ж����ĸ�ά�ȣ�����ж���
	/*printf("termsDimension[%d]:\n", dimensionNo1);
	for (int x = 0; x < termsDimension[dimensionNo1].size(); x++){
		printf("[%s]", termsDimension[dimensionNo1][x].c_str());
	}
	printf("\n");*/
	//û�жԷ�������
	if (literalNo1 == dimensionNo1 * 2){//�����4�������򲻱䣬�������Ӹôʵ�01ֵ
		//�����������1��������Ϊ2�������2��������Ϊ4
		if (termsDimension[dimensionNo1].size() < 4){
			int lsize = termsDimension[dimensionNo1].size();
			for (int x = 0; x < lsize; x++){
				string m;
				if (termsDimension[dimensionNo1][x][0] == '1'){
					m = '0';
				}
				else{
					m = '1';
				}
				m = m + termsDimension[dimensionNo1][x][1];
				termsDimension[dimensionNo1].push_back(m);
				addValues.push_back(m);
			}
		}
	}
	else{
		if (termsDimension[dimensionNo1].size() < 4){
			int lsize = termsDimension[dimensionNo1].size();
			for (int x = 0; x < lsize; x++){
				string m;
				if (termsDimension[dimensionNo1][x][1] == '1'){
					m = '0';
				}
				else{
					m = '1';
				}
				m = termsDimension[dimensionNo1][x][0] + m;
				termsDimension[dimensionNo1].push_back(m);
				addValues.push_back(m);
			}
		}
	}
	//�仯��ĸ񸲸�
	/*printf("\nchanged:\n");
	for (int x = 0; x < dimensionNum; x++){
		for (int y = 0; y < termsDimension[x].size(); y++){
			printf("[%s]", termsDimension[x][y].c_str());
		}
		printf("\n");
	}*/
	//�仯�ĸ񸲸ǵ�
	//printf("\nchanging:\n");
	for (int x = 0; x < dimensionNum; x++){
		if (x == dimensionNo1){
			termsDimension[x].clear();
			termsDimension[x] = addValues;
		}
		/*for (int y = 0; y < termsDimension[x].size(); y++){
			printf("[%s]", termsDimension[x][y].c_str());
		}
		printf("\n");*/
	}
	return termsDimension;
}
bool checkLOFOccur(vector<vector<string>>LOFChangeCells, string model, string testCase){
	for (int i = 0; i < LOFChangeCells.size(); i++){
		string subTestCase = testCase.substr(i * 2, 2);
		//printf("[%s]", subTestCase.c_str());
		bool m = false;
		for (int j = 0; j < LOFChangeCells[i].size(); j++){
			if (subTestCase == LOFChangeCells[i][j]){
				m = true;
				break;
			}
		}
		if (m == false)return false;
	}
	//printf("\n");
	return true;
}


bool checkTestDiffer(string test, string exprOri/*ԭ���ʽ*/, string expr/*����*/){
	Utility uti;
	bool v = (uti.evaluate(exprOri, test)) ^ (uti.evaluateMutant(exprOri, expr, test));//��򣺲�ͬΪ�棬��ͬΪ��
	return v;
}

//��ϡ��������Ų�������
bool diagnosis(
	string oriExp, //ԭ���ʽ
	string faultExp, //������ʽ
	string faultType, //ȱ������
	vector<string> oriTerms,
	vector<string> faultTerms,
	string model,
	vector<vector<string>> optiUniformitySet,//������Ч��
	vector<vector<string>> optiDifferSet//������Ч��
	)
{
	DnfMutant dnfMutant;
	//printf("\nfaultType is [%s]\n", faultType.c_str());
	Mutants z = dnfMutant.getMutants(oriExp, faultType, oriTerms, faultTerms, model);
	if (z.mutants.size() == 0 || z.mutants[0].mutant == ""){
		printf("\nI cannot diagnose the type of faults.\n");
		return false;
	}
	//printf("\nThe number of mutants is [%d]\n", z.mutants.size());
	/*
	for (int i = 0; i < z.mutants.size(); i++){
		printf("%s\n", z.mutants[i].mutant.c_str());
	}
	*/
	

	Mutants realMutants(oriExp);
	for (int j = 0; j < z.mutants.size(); j++){
		string mutant = z.mutants[j].mutant;
		bool flagUniform = true;//��Ч���Ƿ�ȫ��ͨ��
		//��1��������Ч�����̽��
		for (int i = 0; i < optiUniformitySet.size(); i++){
			//printf("point: %s, faultExp: %s, mutant: %s\n", optiUniformitySet[i][0].c_str(), faultExp.c_str(), mutant.c_str());
			if (checkTestDiffer(optiUniformitySet[i][0], faultExp, mutant)){
				flagUniform = false;
				break;
			}
		}

		if (flagUniform == true){//��Ч��ȫ��ͨ��
			bool flagDiffer = true;
			for (int i = 0; i < optiDifferSet.size(); i++){
				if (checkTestDiffer(optiDifferSet[i][0], faultExp, mutant)){
					flagDiffer = false;
					break;
				}
			}
			if (flagDiffer == true){//��Ч��ȫ��ͨ��
				realMutants.addMutant(mutant, z.mutants[j].type, z.mutants[j].oriterms, z.mutants[j].faultterms);
			}
		}
	}
	if (realMutants.mutants.size() == 0){
		return false;
	}
	else{
		//printf("mutant is [%s]   faultExp is [%s]", realMutants.mutants[0].mutant.c_str(), faultType.c_str());
		if (checkEqualment(faultExp, realMutants.mutants[0].mutant)){
			printf("realFaultExpress is[%s][%s][more terms]", realMutants.mutants[0].mutant.c_str(), faultType.c_str());
			return true;
		}
		return true;
	}
}

int main(int argc, char* argv[])
{
	if (argc < 3){
		printf("please enter:cmd expression\nfor example: GlpkDoubleDiagnosis.exe LIORDF ab+cd ab+d");
		return -1;
	}
	string faultType = argv[1];
	string expression = argv[2];
	string faultExpression = "";
	if (argc > 3){
		faultExpression = argv[3];
		printf("faultExpression:%s\n", faultExpression.c_str());
	}
	//�ֽ���ʽ
	PRE_PROCESS_DATA preData;

	clock_t time1 = clock();
	//ת��ΪIDNF
	Utility uti;
	vector<string> terms = uti.extractDNFstyle(expression);
	string model = uti.generateModel(expression);
	string newExp = "";
	//�õ�IDNF
	//newExp = generateIDNForICNF(expression, true);
	//newExp = generateDNF(terms);
	newExp = expression;

	printf("newExp:%s\n", newExp.c_str());

	getPreProcessData(newExp, preData, DNF);
	/*printf("\npreData.op:\n");
	for (int i = 0; i < preData.op.size(); i++){
	printf("%s\n", preData.op[i].c_str());
	}
	printf("\n");
	printf("%d\n", preData.literalCount);
	printf("\npreData.up:\n");
	for (int i = 0; i < preData.up.size(); i++){
	printf("%s\n", preData.up[i].c_str());
	}
	printf("\npreData.termsDimension:\n");
	for (int i = 0; i < preData.termsDimension.size(); i++){
	for (int j = 0; j < preData.termsDimension[i].size(); j++){
	printf("%s\n", preData.termsDimension[i][j].c_str());
	}
	printf("\n");
	}

	printf("\npreData.opover2:\n");
	for (int i = 0; i < preData.opover2.size(); i++){
	printf("%s\n", preData.opover2[i].c_str());
	}
	printf("\n");*/
	CCommonFuncs commonFuncs;
	int literalNum = preData.literalCount;
	int dimensionNum = (literalNum + 1) / 2;
	bool vOdd = (literalNum % 2 == 0) ? false : true;
	//��ȡ���е�OTP
	vector<string> op = preData.op;
	vector<int> allOPs;
	for (int i = 0; i<op.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(op.at(i), vOdd);
		allOPs.push_back(value);
	}
	sort(allOPs.begin(), allOPs.end());

	//��ȡ���е�UTP��
	vector<int> allUPs;
	vector<string> up = preData.up;
	for (int i = 0; i<up.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(up.at(i), vOdd);
		allUPs.push_back(value);
	}
	sort(allUPs.begin(), allUPs.end());

	//��ȡ�Ƿ��е��ʵ��ж�
	bool isSingleLitral = false;
	//����ռȫ��ά����Ŀ
	int fullDemNum = 0;
	//����ռ���ά����Ŀ
	int halfDemNum = 0;
	vector<vector<string>> termsDimension = preData.termsDimension;

	vector<vector<vector<int>>> mutantsDimensionValue;
	vector<vector<CCell*>> corresponding_tests;
	vector<int> tModelIndexAppeared;
	vector<vector<CCell*>> matrixCells;
	vector<vector<int>> coefficientMatrix;
	vector<string> mutants;
	/*printf("%s,%d,termsDimension Print:\n", __FILE__, __LINE__);//�õ���term��ά�ȱ仯��һ��Ϊһ��term�ı仯ά��
	for (int j = 0; j < termsDimension.size(); j++)
	{
	for (int k = 0; k<termsDimension[j].size(); k++)
	printf("%s ", termsDimension[j][k].c_str());
	printf("\n");
	}
	printf("%s,%d,end termsDimension Print\r\n", __FILE__, __LINE__);*/

	vector<vector<CCell *>> TNFCCells;
	vector<string> faultTypes;
	vector<string> faultTerms;

	int startIndex = 0;
	for (int i = 0; i < preData.termsDimension.size(); i++){//���term���н���
		string dimensionValuesStr;
		//��һ��ת��Ϊ��ά�任������11 10 00 01 11 10 ������3 4 1 2 3 4����34 1234
		for (int j = 0; j < dimensionNum; j++)
		{
			//����ά��01��ʾ����ת��Ϊ1,2,3,4��ʾ��
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			int strLength = tempStr.length();
			if (vOdd != true || j != dimensionNum - 1) //����ռȫʱ���ַ�������ӦΪ7,ռ�볤��Ϊ3
			{
				if (strLength == 7)
					fullDemNum++;
				else if (strLength == 3)
					halfDemNum++;
			}
			else //����ռȫʱ���ַ�������ӦΪ3,ռ�볤��Ϊ1
			{
				if (strLength == 3)
					fullDemNum++;
				else if (strLength == 1)
					halfDemNum++;
			}
			if (j != dimensionNum - 1)
				dimensionValuesStr.append("\n");
		}
		if (fullDemNum == dimensionNum - 1 && halfDemNum == 1)
			isSingleLitral = true;
		//����ȱ����Ϣת������Ҫ�ĸ񸲸�,�ڴ˴�Ҫ���и���
		mutantsDimensionValue.clear();
		corresponding_tests.clear();
		//����˼·
		//���ȸ��ݱ���õ����б���Ĵ���㣬��ν����㣬���Ǳ����뱾�岻ͬ����ֵ�㣬������������ֵ�㣬Ҳ����������ֵ��
		//Ȼ����ݴ�����鷽���飬һ������һ�������飬����������ĸ���Ϊ���д����
		//Ȼ�����ÿ������Ĵ����ȷ��ÿ�����̵ı�������
		//�᷽����õ����ŵı����������ɣ����в��漰�����ʽ�����ֵ����������������߼�ԭ����ο����ʽ
		if (faultType == "LIORDF"){
			int matrix_size1 = matrixCells.size();
			LIFGenerateCellsForMutants(commonFuncs,
				dimensionValuesStr,
				vOdd,
				matrixCells,
				mutantsDimensionValue,
				corresponding_tests, allOPs);
			for (int x = startIndex; x < matrixCells.size(); x++){
				faultTypes.push_back("LIF");
				faultTerms.push_back(terms[i]);
			}
			startIndex = matrixCells.size();
			int matrix_size2 = matrixCells.size();
			if (isSingleLitral) //���ֻ��һ���� ���˻���TOF
			{
				TOFGenerateCellsForMutants(commonFuncs, dimensionValuesStr, vOdd, matrixCells,
					mutantsDimensionValue,
					corresponding_tests, allOPs);

				for (int x = startIndex; x < matrixCells.size(); x++){
					faultTypes.push_back("TOF");
					faultTerms.push_back(terms[i]);
				}
				startIndex = matrixCells.size();
			}

			LRFGenerateCellsForMutants(commonFuncs, dimensionValuesStr, vOdd, matrixCells,
				mutantsDimensionValue,
				corresponding_tests, allOPs, allUPs);

			for (int x = startIndex; x < matrixCells.size(); x++){
				faultTypes.push_back("LRF");
				faultTerms.push_back(terms[i]);
			}
			startIndex = matrixCells.size();

			LOFGenerateCellsForMutants(commonFuncs, dimensionValuesStr, vOdd, matrixCells,
				mutantsDimensionValue,
				corresponding_tests, allOPs, allUPs);

			for (int x = startIndex; x < matrixCells.size(); x++){
				faultTypes.push_back("LOF");
				faultTerms.push_back(terms[i]);
			}
			startIndex = matrixCells.size();

			get_mutant_constrains(mutantsDimensionValue, vOdd, corresponding_tests, newExp, i, mutants);

			//���6����CC-S�ɼ���˫ȱ��
			CDoubleLIRFHandler dlirfHandler;
			dlirfHandler.run2(dimensionValuesStr, vOdd, matrixCells, allOPs, allUPs);
			for (int x = startIndex; x < matrixCells.size(); x++){
				mutants.push_back("");
				faultTypes.push_back("LIRF");
				faultTerms.push_back(terms[i]);
			}
			startIndex = matrixCells.size();


			CDoubleLORFHandler dlorfHandler;
			dlorfHandler.run2(dimensionValuesStr, vOdd, matrixCells, allOPs, allUPs);
			for (int x = startIndex; x < matrixCells.size(); x++){
				mutants.push_back("");
				faultTypes.push_back("LORF");
				faultTerms.push_back(terms[i]);
			}
			startIndex = matrixCells.size();


			CDoubleLRRFHandler dlrrfHandler;
			dlrrfHandler.run2(dimensionValuesStr, vOdd, matrixCells, allOPs, allUPs);
			for (int x = startIndex; x < matrixCells.size(); x++){
				mutants.push_back("");
				faultTypes.push_back("LRRF");
				faultTerms.push_back(terms[i]);
			}
			startIndex = matrixCells.size();

		}
		else return 0;

	}

	getLIFxLIFCells2(commonFuncs, terms, dimensionNum, termsDimension, preData.opover2, allOPs, vOdd, matrixCells, mutants, faultTypes, faultTerms);
	getLIFxLRFCells2(commonFuncs, terms, dimensionNum, termsDimension, preData.opover2, allUPs, allOPs, vOdd, matrixCells, mutants, faultTypes, faultTerms);
	getLRFxLRFCells2(commonFuncs, terms, dimensionNum, termsDimension, allUPs, allOPs, vOdd, matrixCells, mutants, faultTypes, faultTerms);

	/*for (int i = 0; i < matrixCells.size(); i++){
	printf("cells[%d]", i);
	for (int j = 0; j < matrixCells[i].size(); j++){
	printf("%d ", matrixCells[i][j]->getTModelIndex());
	}
	printf("\n");
	}*/

	//for (int i = 0; i < mutants.size(); i++){
	//printf("mutants[%d]%s \n", i, mutants[i].c_str());
	//}

	//ת��Լ������
	//�õ����е���ֵ�㣬��һ���������ѭ����
	commonFuncs.getAlltModelIndexAppearedInMatrix(matrixCells, tModelIndexAppeared);
	//������������������311����312֮ǰ
	sort(tModelIndexAppeared.begin(), tModelIndexAppeared.end());

	commonFuncs.getCoefficientMatrix(tModelIndexAppeared, matrixCells, coefficientMatrix);

	vector<vector<int>>finalCoefficientMatrix;
	vector<int>finalTModelIndex;
	removeAllZeroCols(coefficientMatrix, tModelIndexAppeared, finalCoefficientMatrix, finalTModelIndex);

	vector<int>().swap(tModelIndexAppeared);//ֻ��Ϊ�˴�ӡ��
	vector<vector<int>>().swap(coefficientMatrix);//ֻ��Ϊ�˴�ӡ��
	matrixCells.clear();
	//Glpk��Լ������
	GLPK_Solve glpk_Solve;
	vector<int> colSelected;
	vector<int> deletedRows;
	removeAllZeroRows(finalCoefficientMatrix, coefficientMatrix, deletedRows);
	coefficientMatrix.swap(finalCoefficientMatrix);
	coefficientMatrix.clear();
	for (int i = deletedRows.size() - 1; i >= 0; i--){
		vector<string>::iterator itrl = mutants.begin() + i;
		mutants.erase(itrl);
	}
	//�����޸ģ�����Ҫÿ�ζ�����
	vector<string> mutants_DP;
	for (int num = 0; num<mutants.size(); num++)
	{
		//printf("%s,%d,mutants[%d]:[%s]\r\n", __FILE__, __LINE__, num, mutants[num].c_str());
		mutants_DP.push_back(mutants[num]);
	}

	glpk_Solve.run(colSelected, &finalCoefficientMatrix);
	//���õ��ĸ񸲸�ת���ɲ�������
	vector<string> testset;
	for (int j = 0; j < colSelected.size(); j++)
	{
		char tmp[40];
		sprintf(tmp, "%d", finalTModelIndex[colSelected[j]]);
		//printf("%d:", finalTModelIndex[colSelected[j]]);
		string inde = "";
		inde.append(tmp);
		string point01Str = commonFuncs.transTModelIndexStrTo01Str(inde, vOdd);
		testset.push_back(point01Str);
		//printf("%s\n", point01Str.c_str());
	}
	clock_t time2 = clock();
	//printf("testcase num[%d] time[%d] mutants num[%d]\n", testset.size(), time2 - time1, mutants.size());
	//printf("���������������[%s]\n", faultExpression.c_str());

	vector<vector<int>> DP_Matrix;
	vector<vector<int>> case_mutant;

	genDP_Matrix(finalCoefficientMatrix, finalTModelIndex, DP_Matrix);//��Լ�����������ת������һ�оͶ�Ӧһ���񸲸�

	genCase_Mutants(case_mutant, finalTModelIndex, DP_Matrix);
	vector<int> finalTModelIndexMatrixCPY(finalTModelIndex);
	vector<vector<int>> finalCoefficientMatrixCPY(finalCoefficientMatrix);
	FaultDetection faultDetector;
	if (faultExpression == ""){
		//�����ʽ
		Mutants zMutants(expression);
		DnfMutant dnfMutant;
		zMutants = dnfMutant.getMutants(newExp, faultType);

		for (int l = 0; l < zMutants.mutants.size(); l++)
		{
			string result2;
			result2 = faultDetector.checkExprwithSet(newExp,
				zMutants.mutants[l].mutant,
				testset,
				case_mutant,
				finalCoefficientMatrix,
				finalTModelIndex,
				mutants_DP,
				vOdd);
			if (result2 != "no fault")
			{
				if (faultDetector.checkEqualment(result2, zMutants.mutants[l].mutant))
				{
					//printf("enter checkEqualment true\n");
				}
				else
				{
					//printf("enter checkEqualment false\n");
				}
			}
			else
			{
				//printf("no fault\n");
			}
		}
	}
	else{
		printf("OriExp:%s, faultExp:%s\n", newExp.c_str(), faultExpression.c_str());
		string result2;
		vector<vector<string>>maybeMutants;
		vector<vector<string>> uniformitySet;
		vector<vector<string>> differSet;
		vector<vector<string>> optiUniformitySet;
		vector<vector<string>> optiDifferSet;
		int  isFirstO = true;

		vector<int> selectAddRows;
		{//����ѭ���ң�ֱ���ҷ��̼���
			vector<int> uniformityTCsol, differTCsol;
			vector<string> uniformityStr, differStr;
			vector<string> uniformityType, differType;
			//��ȡ��Ч�����Ч��
			if (checkTestSet(
				testset,
				newExp, faultExpression,
				uniformityTCsol, differTCsol, uniformityStr, differStr, uniformityType, differType, vOdd))
			{//�����棬˵��ȫ����Ч���Ե�
				printf("\ncannot find any faults\n");
				return 0;
			}
			else{
				bool isExpand = false;
				bool isShrink = false;
				//��¼��Ч��
				for (int l = 0; l < uniformityStr.size(); l++){
					vector<string> llm;
					llm.push_back(uniformityStr[l]);
					llm.push_back(uniformityType[l]);
					uniformitySet.push_back(llm);
					if (isFirstO){
						optiUniformitySet.push_back(llm);
					}
				}
				//��¼��Ч��
				for (int l = 0; l < differStr.size(); l++){
					if (differType[l] == "e"){//����������
						isExpand = true;
					}
					else{
						isShrink = true;
					}
					vector<string> llm;
					llm.push_back(differStr[l]);
					llm.push_back(differType[l]);
					differSet.push_back(llm);
					if (isFirstO){
						optiDifferSet.push_back(llm);
					}
				}
				isFirstO = false;

				if (uniformityTCsol.size() != 0)sort(uniformityTCsol.begin(), uniformityTCsol.end());//һ�µĵ�
				if (differTCsol.size() != 0)sort(differTCsol.begin(), differTCsol.end());//��һ�µĵ�
				/*printf("The UniformityTCsol is \n");
				for (int i = 0; i < uniformityTCsol.size(); i++){
					printf("[%d]", uniformityTCsol[i]);
				}
				printf("\n");

				printf("The DifferTCsol is \n");
				for (int i = 0; i < differTCsol.size(); i++){
					printf("[%d]", differTCsol[i]);
				}
				printf("\n");

				printf("finalCoefficientMatrix:\n");
				for (int i = 0; i < finalTModelIndex.size(); i++){
					printf("[%d]", finalTModelIndex[i]);
				}
				printf("\n");
				for (int i = 0; i < finalCoefficientMatrix.size(); i++){
					for (int j = 0; j < finalCoefficientMatrix[i].size(); j++){
						printf("[%d] ", finalCoefficientMatrix[i][j]);
					}
					printf("\n");
				}
				printf("\n");


				for (int j = 0; j < finalCoefficientMatrix.size(); j++)
				{
					printf("mutant[%d][%s][%s][%s]\n", j, mutants[j].c_str(), faultTypes[j].c_str(), faultTerms[j].c_str());
				}
				printf("\n");

				printf("DP_Matrix:\n");
				for (int i = 0; i < DP_Matrix.size(); i++){
					for (int j = 0; j < DP_Matrix[i].size(); j++){
						printf("[%d]", DP_Matrix[i][j]);
					}
					printf("\n");
				}
				printf("\n");*/
				vector<int> colSelected1;
				glpk_Solve.run_lo(colSelected1, &DP_Matrix);//������ʵ����ɾ������������ת��⣬�����һ���ģ�����ʦû˵��������û�����۾�û���˺�
				/*printf("colSelected1 1211:\n");
				for (int i = 0; i < colSelected1.size(); i++){
					printf("[%d]", colSelected1[i]);
				}
				printf("\n");*/
				for (int j = 0; j < colSelected1.size(); j++)
				{
					//printf("mutant[%d][%s][%s]\n", colSelected1[j], mutants[colSelected1[j]].c_str(), faultTypes[colSelected1[j]].c_str());
					//�����LIF�������⴦��ֻ���������򣬸ñ��ʽ��������������Ӧ��LRF����������ŵ㣬�򲻱����ñ��ʽ��������Ӧ��LRF��ֻ�����ŵ㣬��ֱ�ӿ�LOF
					if (faultTypes[colSelected1[j]] == "LIF"){
						bool remainSign = false;
						bool tmt = true;
						if (isExpand){//ֻ�����žͲ��Ƚ��ˣ�ֱ����һ��
							tmt = false;
						}
						else{
							for (int y = 0; y < uniformityStr.size(); y++){//��Ч����������
								bool v2 = (uti.evaluateMutant(newExp, mutants[colSelected1[j]], uniformityStr[y]));
								if (uti.evaluate(newExp, uniformityStr[y]) ^ v2){//����������
									tmt = false;
									break;
								}
							}
						}
						if (tmt){
							//printf("LIF remain\n");
							remainSign = true;//���е㶼��ȷ������
						}
						if (remainSign){
							//printf("maybe a mutant[%d][%s][%s][%s]\n", colSelected1[j], mutants[colSelected1[j]].c_str(), faultTypes[colSelected1[j]].c_str(), faultTerms[colSelected1[j]].c_str());
							vector<string> ltsd;
							ltsd.push_back(mutants[colSelected1[j]]);
							ltsd.push_back(faultTypes[colSelected1[j]]);
							ltsd.push_back(faultTerms[colSelected1[j]]);
							maybeMutants.push_back(ltsd);
							selectAddRows.push_back(colSelected1[j]);
						}
						{//��LRF
							//ֱ���Ҷ�Ӧ��LRF�����ɷ�����ֱ��LOF�����޳��ȼ۱��ʽ
							DnfMutant dnfMutant;
							vector<string> faultExpTerms = uti.extractDNFstyle(mutants[colSelected1[j]]);
							vector<string> faultOccurTerms;
							set_difference(
								terms.begin(),
								terms.end(),
								faultExpTerms.begin(),
								faultExpTerms.end(),
								inserter(
								faultOccurTerms,
								faultOccurTerms.end()
								)
								);
							for (int x = 0; x < faultOccurTerms.size(); x++){
								if (faultOccurTerms[x].size() == 1 || faultOccurTerms[x].size() == 2 && faultOccurTerms[x][0] == '!'){
									vector<string>::iterator itr = faultOccurTerms.begin() + x;
									faultOccurTerms.erase(itr);
								}
							}
							if (faultOccurTerms.size() == 0)continue;
							Mutants ab = dnfMutant.generateLOFmutants(mutants[colSelected1[j]], faultOccurTerms);//LIF&LOF����LRF
							for (int y = 0; y < ab.mutants.size(); y++){
								//����ñ��ʽ��ԭ���ʽ�ȼۣ����޳�
								if (checkEqualment(newExp, ab.mutants[y].mutant)){
									continue;
								}
								//printf("Find LRF mutant[%s][%s][%s]\n", ab.mutants[y].mutant.c_str(), "LRF", faultTerms[colSelected1[j]].c_str());
								//������������������Ƿ���ȷ
								bool tmt = true;
								for (int z = 0; z < uniformityStr.size(); z++){
									bool v2 = (uti.evaluateMutant(newExp, ab.mutants[y].mutant, uniformityStr[z]));
									if (uti.evaluate(newExp, uniformityStr[z]) ^ v2){//����������
										tmt = false;
										break;
									}
								}
								if (tmt){//���е㶼��ȷ������
									//printf("maybe a LRF mutant[%s][%s][%s]\n", ab.mutants[y].mutant.c_str(), "LRF", faultTerms[colSelected1[j]].c_str());
									vector<string> ltsd;
									ltsd.push_back(ab.mutants[y].mutant);
									ltsd.push_back("LRF");
									ltsd.push_back(faultTerms[colSelected1[j]]);
									maybeMutants.push_back(ltsd);
								}
							}
						}
					}
					else{
						vector<int> uniformityTCsol1, differTCsol1;
						vector<string> uniformityStr1, differStr1;
						vector<string> uniformityType1, differType1;

						if (mutants[colSelected1[j]] == ""){//û�����ɱ��壬��ֱ�Ӽ���Լ������
							/*���ٿ����ˣ�ֱ���õ�ȱ���Ʋ�˫ȱ��*/
							int uvl = finalCoefficientMatrix[colSelected1[j]].size();
							int uvn = finalTModelIndex.size();
							bool rightSign = true;
							for (int x = 0; x < uniformityTCsol.size(); x++){
								for (int y = 0; y < finalTModelIndex.size(); y++){
									if (uniformityTCsol[x] == finalTModelIndex[y]){
										if (finalCoefficientMatrix[colSelected1[j]][y] != 0){//����������
											rightSign = false;
											break;
										}
									}
								}
								if (!rightSign)break;
							}
							if (!rightSign)continue;
							else{//�ҵ���һ��Լ������
								//printf("maybe a mutant[%d][%s][%s][%s]\n", colSelected1[j], mutants[colSelected1[j]].c_str(), faultTypes[colSelected1[j]].c_str(), faultTerms[colSelected1[j]].c_str());
								if (find(selectAddRows.begin(), selectAddRows.end(), colSelected1[j]) == selectAddRows.end()){
									vector<string> ltsd;
									char lz[10];
									sprintf(lz, "#%d", colSelected1[j]);
									ltsd.push_back(lz);//ֱ�Ӽ�¼��ź���
									ltsd.push_back(faultTypes[colSelected1[j]]);
									ltsd.push_back(faultTerms[colSelected1[j]]);
									maybeMutants.push_back(ltsd);
									selectAddRows.push_back(colSelected1[j]);

								}
							}
						}
						else{
							if (checkTestSet(
								testset,
								newExp, mutants[colSelected1[j]],
								uniformityTCsol1, differTCsol1, uniformityStr1, differStr1, uniformityType1, differType1, vOdd))
							{//�����棬˵��ȫ����Ч���Ե�
								continue;
							}
							else{
								bool rightSign = true;
								if (faultTypes[colSelected1[j]] == "LOF"){
									//����������㣬ֱ���ų�
									if (isShrink){
										rightSign = false;
									}
									else{
										for (int x = 0; x < uniformityTCsol.size(); x++){
											bool v2 = (uti.evaluateMutant(newExp, mutants[colSelected1[j]], uniformityStr[x]));
											if (uti.evaluate(newExp, uniformityStr[x]) ^ v2){//����������
												rightSign = false;
												break;
											}
										}
									}
								}
								else{
									for (int x = 0; x < uniformityTCsol.size(); x++){
										bool v2 = (uti.evaluateMutant(newExp, mutants[colSelected1[j]], uniformityStr[x]));
										if (uti.evaluate(newExp, uniformityStr[x]) ^ v2){//����������
											rightSign = false;
											break;
										}
									}
								}
								if (!rightSign)continue;
								else{//�ҵ���һ��Լ������
									//printf("maybe a mutant[%d][%s][%s][%s]\n", colSelected1[j], mutants[colSelected1[j]].c_str(), faultTypes[colSelected1[j]].c_str(), faultTerms[colSelected1[j]].c_str());
									vector<string> ltsd;
									ltsd.push_back(mutants[colSelected1[j]]);
									ltsd.push_back(faultTypes[colSelected1[j]]);
									ltsd.push_back(faultTerms[colSelected1[j]]);
									maybeMutants.push_back(ltsd);
									selectAddRows.push_back(colSelected1[j]);
								}
							}
						}
					}
				}
				//����������������LOF���������͵ı���򷽳̶���¼�ڰ���ʹ����Ч���ж�
				for (int x = 0; x < finalCoefficientMatrix.size(); x++){
					if (find(selectAddRows.begin(), selectAddRows.end(), x) == selectAddRows.end()){//˵��û�иü�¼
						bool ntl = true;
						for (int i = 0; i < uniformityTCsol.size() && ntl; i++){
							for (int j = 0; j < finalTModelIndex.size() && ntl; j++){
								if (uniformityTCsol[i] == finalTModelIndex[j]){//�ҵ������������������Щ���̿��Թҹ�
									if (finalCoefficientMatrix[x][j] != 0){//��Ч��
										ntl = false;
										break;
									}
								}
							}
						}
						if (ntl){
							vector<string> ltsd;
							if (mutants[x] != "")ltsd.push_back(mutants[x]);
							else{
								/*char lz[10];
								sprintf(lz, "#%d", x);
								ltsd.push_back(lz);*/
								continue;
							}
							ltsd.push_back(faultTypes[x]);
							ltsd.push_back(faultTerms[x]);
							maybeMutants.push_back(ltsd);
							selectAddRows.push_back(x);
						}
					}
				}
			}
		}
		//printf("\n");
		//����Ƿ���û����֤�Ĳ�������������У���ֱ�ӽ��м��飬������֤����
		/*printf("finalTModelIndex 1213:\n");
		for (int j = 0; j < finalTModelIndex.size(); j++)
		{
			printf("[%d]", finalTModelIndex[j]);
		}
		printf("\n");*/
		testset.clear();
		for (int j = 0; j < finalTModelIndex.size(); j++)
		{
			char tmp[40];
			sprintf(tmp, "%d", finalTModelIndex[j]);
			//printf("%d:", finalTModelIndex[j]);
			string inde = "";
			inde.append(tmp);
			string point01Str = commonFuncs.transTModelIndexStrTo01Str(inde, vOdd);
			testset.push_back(point01Str);
			//printf("%s\n", point01Str.c_str());
		}
		vector<int> uniformityTCsol, differTCsol;
		vector<string> uniformityStr, differStr;
		vector<string> uniformityType, differType;
		//��ȡ��Ч�����Ч��
		checkTestSet(
			testset,
			newExp, faultExpression,
			uniformityTCsol, differTCsol, uniformityStr, differStr, uniformityType, differType, vOdd);

		//������Ը������������ж����ŵ��Ӧ��term
		vector<string> shrinkTestCase;
		for (int i = 0; i < differType.size(); i++){
			if (differType[i] == "s")shrinkTestCase.push_back(differStr[i]);
		}
		vector<int> shrinkTermsNo = checkShrinkTerms(newExp, terms, shrinkTestCase);

		/*
		for (int i = 0; i < shrinkTermsNo.size(); i++){
			printf("shrinkTerms[%d][%s]\n", shrinkTermsNo[i], terms[shrinkTermsNo[i]].c_str());
		}
		*/

		vector<string> expandTestCase;
		for (int i = 0; i < differType.size(); i++){
			if (differType[i] == "e")expandTestCase.push_back(differStr[i]);
		}
		//����ȷ���Ƿ������ŵ㣬����У�ȷ��������������ŵ�
		vector<vector<string>> expandTermsByETerms;
		for (int i = 0; i < maybeMutants.size(); i++){
			if (maybeMutants[i][1] == "LOF"){
				for (int j = 0; j < maybeMutants[i][2].size(); j++){
					int expandTermsNo = getTermNo(terms, maybeMutants[i][2]);
					string Liter1;
					if (maybeMutants[i][2][j] != '!')Liter1 = maybeMutants[i][2][j];
					else{
						Liter1 = "!";
						j++;
						Liter1 += maybeMutants[i][2][j];
					}
					vector<vector<string>> LOFChangeCells = getLOFChangeCells(expression, model, expandTermsNo, Liter1);
					for (int k = expandTestCase.size() - 1; k >= 0; k--){
						//printf("[%s]:\n", expandTestCase[k].c_str());
						bool nlt = checkLOFOccur(LOFChangeCells, model, expandTestCase[k]);
						if (nlt){
							//printf("belong the expand term[%d]\n", expandTermsNo);
							vector<string>::iterator itr = expandTestCase.begin() + k;
							expandTestCase.erase(itr);
						}
						else{
							//printf("belong the expand term[%d]\n", expandTermsNo);
							continue;
						}
					}
					LOFChangeCells.clear();
				}
			}
		}
		/*printf("remain expand case:");
		for (int i = 0; i < expandTestCase.size(); i++){
			printf("[%s]", expandTestCase[i].c_str());
		}
		printf("\n");*/

		vector<vector<string>> expandTermsByTestCase;
		vector<string>shrinkTerms;
		for (int i = 0; i< shrinkTermsNo.size(); i++){
			shrinkTerms.push_back(terms[shrinkTermsNo[i]]);
		}

		for (int i = shrinkTermsNo.size() - 1; i >= 0; i--){//�ж����������Ƿ������ŵ�
			bool lnd = false;
			for (int j = 0; j < terms[shrinkTermsNo[i]].size(); j++){
				string Liter1;
				if (terms[shrinkTermsNo[i]][j] != '!')Liter1 = terms[shrinkTermsNo[i]][j];
				else{
					Liter1 = "!";
					j++;
					Liter1 += terms[shrinkTermsNo[i]][j];
				}
				vector<vector<string>> LOFChangeCells = getLOFChangeCells(expression, model, shrinkTermsNo[i], Liter1);
				for (int k = 0; k < expandTestCase.size(); k++){
					//printf("[%s]:\n", expandTestCase[k].c_str());
					bool nlt = checkLOFOccur(LOFChangeCells, model, expandTestCase[k]);
					if (nlt){
						//printf("belong the term[%d]\n", shrinkTermsNo[i]);
						if (expandTermsByTestCase.size() <= k){
							vector<string> a;
							a.push_back(terms[shrinkTermsNo[i]]);
							expandTermsByTestCase.push_back(a);
							a.clear();
							lnd = true;
						}
						else{
							expandTermsByTestCase[k].push_back(terms[shrinkTermsNo[i]]);
							lnd = true;
						}
					}
					else{
						//printf("not belong the term[%d]\n", shrinkTermsNo[i]);
						if (expandTermsByTestCase.size() <= k){
							vector<string> a;
							expandTermsByTestCase.push_back(a);
						}
						continue;
					}
				}
				/*for (int k = 0; k < expandTermsByTestCase.size(); k++){
					printf("[%s]:", expandTestCase[k].c_str());
					for (int l = 0; l < expandTermsByTestCase[k].size(); l++){
						printf("[%s]", expandTermsByTestCase[k][l].c_str());
					}
					printf("\n");
				}
				printf("\n");*/
				LOFChangeCells.clear();
			}
			/*if (lnd){
				vector<string>::iterator itrl = shrinkTerms.begin() + i;
				shrinkTerms.erase(itrl);
			}*/
		}
		vector<string> shrinkExpandTerms;
		for (int i = 0; i < expandTermsByTestCase.size(); i++){
			for (int j = 0; j < expandTermsByTestCase[i].size(); j++){
				if (expandTermsByTestCase[i][j] != ""){
					if (find(shrinkExpandTerms.begin(), shrinkExpandTerms.end(), expandTermsByTestCase[i][j]) == shrinkExpandTerms.end())
						shrinkExpandTerms.push_back(expandTermsByTestCase[i][j]);
				}
			}
		}

		/*
		for (int i = 0; i < shrinkExpandTerms.size(); i++){
			printf("shrinkExpandTerms[%d]:%s\n", i, shrinkExpandTerms[i].c_str());
		}
		*/

		//��¼��Ч��
		for (int l = 0; l < uniformityStr.size(); l++){
			vector<string> llm;
			llm.push_back(uniformityStr[l]);
			llm.push_back(uniformityType[l]);
			uniformitySet.push_back(llm);
		}
		//��¼��Ч��
		for (int l = 0; l < differStr.size(); l++){
			vector<string> llm;
			llm.push_back(differStr[l]);
			llm.push_back(differType[l]);
			differSet.push_back(llm);
		}

		/*
		printf("maybeMutants 1213:\n");
		for (int i = 0; i < maybeMutants.size(); i++){
			printf("[%s][%s][%s]\n", maybeMutants[i][0].c_str(), maybeMutants[i][1].c_str(), maybeMutants[i][2].c_str());
		}
		printf("\n");
		*/

		//ʹ����Ч����������Ź���maybeMutants
		bool sSign = false;
		bool eSign = false;
		for (int j = 0; j < differSet.size(); j++){
			if (differSet[j][1] == "s")sSign = true;
			if (differSet[j][1] == "e")eSign = true;
			if (sSign && eSign)break;
		}
		for (int i = maybeMutants.size() - 1; i >= 0;){
			bool lmt = false;
			if (sSign && !eSign){//ֻ���������޳�LOF
				if (maybeMutants[i][1] == "LOF" || maybeMutants[i][1] == "LRF"){//��ʱLRF��LIF�ĵȼ۱��ʽ
					vector<vector<string>>::iterator itr = maybeMutants.begin() + i;
					maybeMutants.erase(itr);
					i--;
					if (i < 0)break;
					lmt = true;
				}
			}
			if (eSign && !sSign){//ֻ�����ţ��޳�LIF
				if (maybeMutants[i][1] == "LIF" || maybeMutants[i][1] == "TOF"){
					vector<vector<string>>::iterator itr = maybeMutants.begin() + i;
					maybeMutants.erase(itr);
					i--;
					if (i < 0)break;
					lmt = true;
				}
			}
			if (!lmt)i--;
		}
		//ʹ����Ч����ȷ�Թ���maybeMutants
		for (int i = maybeMutants.size() - 1; i >= 0; i--){
			bool selectedSign = true;//���߼��������˫ȱ�ݣ����ٿ���û�б���ķ���
			if (maybeMutants[i][0][0] == '#'){//û�б���
				string l = maybeMutants[i][0].substr(1, maybeMutants[i][0].size() - 1);
				int ln = atoi(l.c_str());
				for (int i = 0; i < uniformityTCsol.size() && selectedSign; i++){
					for (int j = 0; j < finalTModelIndex.size() && selectedSign; j++){
						if (uniformityTCsol[i] == finalTModelIndex[j]){//�ҵ������������������Щ���̿��Թҹ�
							if (finalCoefficientMatrix[ln][j] != 0){//��Ч��
								selectedSign = false;
								break;
							}
						}
					}
				}
			}
			else{
				for (int j = 0; j < uniformitySet.size(); j++){
					selectedSign = true;
					bool v1 = (uti.evaluate(newExp, uniformitySet[j][0]));
					bool v2 = (uti.evaluateMutant(newExp, maybeMutants[i][0], uniformitySet[j][0]));
					if (v1 ^ v2){//ֻҪ��һ�����Ծ��޳�
						selectedSign = false;
						break;
					}
				}
				//ʹ�ñ������Ч�����
				/*vector<int> uniformityTCsol1, differTCsol1;
				vector<string> uniformityStr1, differStr1;
				vector<string> uniformityType1, differType1;
				checkTestSet(
				testset,
				newExp, maybeMutants[i][0],
				uniformityTCsol1, differTCsol1, uniformityStr1, differStr1, uniformityType1, differType1, vOdd);
				for (int j = 0; j < differStr1.size(); j++){
				selectedSign = true;
				bool v1 = (uti.evaluate(newExp, differStr1[j]));
				bool v2 = (uti.evaluateMutant(newExp, maybeMutants[i][0], differStr1[j]));
				if (!(v1 ^ v2)){//ֻҪ��һ�����Ծ��޳�
				selectedSign = false;
				break;
				}
				}*/
			}
			if (!selectedSign){
				vector<vector<string>>::iterator itr = maybeMutants.begin() + i;
				maybeMutants.erase(itr);
			}
		}

		/*
		printf("realEntranceMutants 1214:\n");
		for (int i = 0; i < maybeMutants.size(); i++){
			printf("[%s][%s][%s]\n", maybeMutants[i][0].c_str(), maybeMutants[i][1].c_str(), maybeMutants[i][2].c_str());
		}
		printf("\n");
		*/

		for (int i = maybeMutants.size() - 1; i >= 0; i--){
			if (maybeMutants[i][1] == "LIF"){
				vector<vector<string>>::iterator itrl = maybeMutants.begin() + i;
				maybeMutants.erase(itrl);
			}
		}

		//���maybeMutants��û�����ݣ�˵���϶���˫ȱ��
		//�����������ΪLIF��LOF������term����ʾ�϶��ǲ���LRF��˫��˫ȱ��
		//�����������ΪLIF��LOF��һ��term����ʾ�϶�Ϊ��ȱ�ݣ��򲻺�LRF�ĵ���˫ȱ��
		//�����������ΪLRF��һ��term����ʾ�϶�Ϊ��ȱ�ݣ���LRF�ĵ���˫ȱ��
		//���û�����ݣ�����÷�����Ч�����,����Ч�����
		if (maybeMutants.size() == 0){
			for (int i = maybeMutants.size() - 1; i >= 0; i--){
				bool selectedSign = true;
				for (int j = 0; j < differSet.size(); j++){
					bool v1 = (uti.evaluate(newExp, differSet[j][0]));
					bool v2 = (uti.evaluateMutant(newExp, maybeMutants[i][0], differSet[j][0]));
					if (v1 ^ v2){//ֻҪ��һ�����Ծ��޳�
						selectedSign = false;
						break;
					}
				}
				if (!selectedSign){
					vector<vector<string>>::iterator itr = maybeMutants.begin() + i;
					maybeMutants.erase(itr);
				}
			}

			/*
			printf("realEntranceMutants 1215:\n");
			for (int i = 0; i < maybeMutants.size(); i++){
				printf("[%s][%s][%s]\n", maybeMutants[i][0].c_str(), maybeMutants[i][1].c_str(), maybeMutants[i][2].c_str());
			}
			printf("\n");
			*/

		}

		//������ʽ��������Ҫ�ı��ʽ����ֱ���˳���
		for (int i = 0; i < maybeMutants.size(); i++){
			if (maybeMutants[i][0][0] == '#'){
			}
			else{
				if (checkEqualment(faultExpression, maybeMutants[i][0])){
					printf("realFaultExpress is[%s][%s][%s]\n", maybeMutants[i][0].c_str(), maybeMutants[i][1].c_str(), maybeMutants[i][2].c_str());
					return 0;
				}
			}
		}
		//��ȱ�����ͺͷ���ȱ������й���
		vector<string>maybeTypes;
		vector<vector<string>>maybeTerms;
		//�鲢ȱ������
		for (int i = 0; i < maybeMutants.size(); i++){
			bool lmdt = true;
			for (int j = 0; j < maybeTypes.size(); j++){
				if (maybeMutants[i][1] == maybeTypes[j]){//�ҵ���ȱ���������ò���
					lmdt = false;
					break;
				}
			}
			if (lmdt){
				maybeTypes.push_back(maybeMutants[i][1]);
				vector<string> mm;
				maybeTerms.push_back(mm);
			}
		}
		//�鲢ȱ�ݷ�����
		for (int i = 0; i < maybeMutants.size(); i++){
			bool lmdt = true;
			for (int j = 0; j < maybeTypes.size(); j++){
				if (maybeMutants[i][1] == maybeTypes[j]){//�ҵ���ȱ������
					bool muv = true;
					for (int k = 0; k < maybeTerms[j].size(); k++){
						if (maybeTerms[j][k] == maybeMutants[i][2]){//�����Ѿ�����
							muv = false;
						}
					}
					if (muv){//�������
						maybeTerms[j].push_back(maybeMutants[i][2]);
					}
					break;
				}
			}
		}

		bool nldt = false;
		for (int i = 0; i < maybeTypes.size(); i++){
			if (maybeTypes[i] == "LRF"){
				nldt = true;
				for (int j = 0; j < shrinkExpandTerms.size(); j++){
					if (find(maybeTerms[i].begin(), maybeTerms[i].end(), shrinkExpandTerms[j]) == maybeTerms[i].end()){
						maybeTerms[i].push_back(shrinkExpandTerms[j]);
					}
				}
				break;
			}
		}
		if (!nldt && shrinkExpandTerms.size() > 0){
			maybeTypes.push_back("LRF");
			vector<string> ll;
			for (int j = 0; j < shrinkExpandTerms.size(); j++){
				ll.push_back(shrinkExpandTerms[j]);
			}
			maybeTerms.push_back(ll);
		}

		nldt = false;
		for (int i = 0; i < maybeTypes.size(); i++){
			if (maybeTypes[i] == "LIF"){
				nldt = true;
				for (int j = 0; j < shrinkTerms.size(); j++){
					if (find(maybeTerms[i].begin(), maybeTerms[i].end(), shrinkTerms[j]) == maybeTerms[i].end()){
						maybeTerms[i].push_back(shrinkTerms[j]);
					}
				}
				break;
			}
		}
		if (!nldt && shrinkTerms.size() > 0){
			maybeTypes.push_back("LIF");
			vector<string> ll;
			for (int j = 0; j < shrinkTerms.size(); j++){
				ll.push_back(shrinkTerms[j]);
			}
			maybeTerms.push_back(ll);
		}

		//���������û��LIF��ֱ���޳�LRF�����ֻ��һ��LIF��ֱ���޳�LRFxLRF���ж�����е�LRF�е�������
		/*
		printf("types and terms 1212:[%d]\n", maybeTypes.size());
		for (int i = 0; i < maybeTypes.size(); i++){
			printf("[%s]\n", maybeTypes[i].c_str());
			for (int j = 0; j < maybeTerms[i].size(); j++){
				printf("[%s]", maybeTerms[i][j].c_str());
			}
			printf("\n\n");
		}
		printf("\n\n");
		*/


		vector<string> lifTerms;
		bool existLif = false;
		for (int i = 0; i < maybeTypes.size(); i++){
			if (maybeTypes[i] == "LIF"){
				existLif = true;
				lifTerms = maybeTerms[i];
			}
		}
		if (existLif == false){
			//�޳�LRF��LRFxLRF
			for (int i = 0; i < maybeTypes.size(); i++){
				if (maybeTypes[i] == "LRF"){
					vector<string>::iterator itrn = maybeTypes.begin() + i;
					maybeTypes.erase(itrn);
					vector<vector<string>>::iterator itrm = maybeTerms.begin() + i;
					maybeTerms.erase(itrm);
				}
				if (maybeTypes[i] == "LRFxLRF"){
					vector<string>::iterator itrn = maybeTypes.begin() + i;
					maybeTypes.erase(itrn);
					vector<vector<string>>::iterator itrm = maybeTerms.begin() + i;
					maybeTerms.erase(itrm);
				}
			}
		}
		else{
			if (lifTerms.size() == 1){
				//�޳�LRFxLRF��LRF�����е�������
				for (int i = 0; i < maybeTypes.size(); i++){
					if (maybeTypes[i] == "LRF"){
						for (int j = 0; j < maybeTerms[i].size(); j++){
							if (maybeTerms[i][j] != lifTerms[0]){
								vector<string>::iterator itrn = maybeTerms[i].begin() + j;
								maybeTerms[i].erase(itrn);
							}
						}
						if (maybeTerms[i].size() == 0){
							vector<string>::iterator itrn = maybeTypes.begin() + i;
							maybeTypes.erase(itrn);
							vector<vector<string>>::iterator itrm = maybeTerms.begin() + i;
							maybeTerms.erase(itrm);
						}
					}
					if (maybeTypes[i] == "LRFxLRF"){
						vector<string>::iterator itrn = maybeTypes.begin() + i;
						maybeTypes.erase(itrn);
						vector<vector<string>>::iterator itrm = maybeTerms.begin() + i;
						maybeTerms.erase(itrm);
					}
				}
			}
		}

		for (int i = 0; i < maybeTypes.size(); i++){
			for (int j = 0; j < maybeTerms[i].size(); j++){
				string::iterator itrl = find(maybeTerms[i][j].begin(), maybeTerms[i][j].end(), '+');
				if (itrl != maybeTerms[i][j].end()){
					vector<string> bTerms = uti.extractDNFstyle(maybeTerms[i][j]);
					maybeTerms[i][j] = bTerms[0];
					for (int k = 1; k < bTerms.size(); k++){
						maybeTerms[i].push_back(bTerms[k]);
					}
				}
			}
		}

		/*
		printf("types and terms 1213:[%d]\n", maybeTypes.size());
		for (int i = 0; i < maybeTypes.size(); i++){
			printf("[%s]\n", maybeTypes[i].c_str());
			for (int j = 0; j < maybeTerms[i].size(); j++){
				printf("[%s]", maybeTerms[i][j].c_str());
			}
			printf("\n\n");
		}
		printf("\n\n");
		*/

		vector<vector<string>> termCharacters;
		vector<string> termsChanged;
		for (int i = 0; i < maybeTerms.size(); i++){
			for (int j = 0; j < maybeTerms[i].size(); j++){
				vector<string>::iterator itr = find(termsChanged.begin(), termsChanged.end(), maybeTerms[i][j]);
				if (itr == termsChanged.end()){
					termsChanged.push_back(maybeTerms[i][j]);
					vector<string> a;
					a.push_back(maybeTerms[i][j]);
					if (maybeTypes[i] == "LIF" || maybeTypes[i] == "TOF" || maybeTypes[i] == "LIFxLIF"){
						a.push_back("shrink");
					}
					else if (maybeTypes[i] == "LOF"){
						a.push_back("expand");
					}
					else if (maybeTypes[i] == "LRF" || maybeTypes[i] == "LRFxLRF"){
						a.push_back("shrink+expand");
					}
					else{
						a.push_back("shrink+expand");
					}
					termCharacters.push_back(a);
				}
				else{
					int ln = itr - termsChanged.begin();
					if (maybeTypes[i] == "LIF" || maybeTypes[i] == "TOF" || maybeTypes[i] == "LIFxLIF"){
						if (termCharacters[ln][1] == "expand")termCharacters[ln][1] = "shrink+expand";
					}
					else if (maybeTypes[i] == "LOF"){
						if (termCharacters[ln][1] == "shrink")termCharacters[ln][1] = "shrink+expand";
					}
					else if (maybeTypes[i] == "LRF" || maybeTypes[i] == "LRFxLRF"){
						termCharacters[ln][1] = "shrink+expand";
					}
					else{
						termCharacters[ln][1] = "shrink+expand";
					}
				}
			}
		}

		//�����ӿڴ�ӡ��
		printf("terms and types 1213-1:[%d]\n", termCharacters.size());
		for (int i = 0; i < termCharacters.size(); i++){
			for (int j = 0; j < termCharacters[i].size(); j++){
				printf("[%s]", termCharacters[i][j].c_str());
			}
			printf("\n");
		}
		printf("\n\n");

		/*
		printf("differSet 1213:\n");
		for (int i = 0; i < differSet.size(); i++){
			printf("[%s][%s]\n", differSet[i][0].c_str(), differSet[i][1].c_str());
		}
		printf("\n");

		printf("uniformitySet 1213:\n");
		for (int i = 0; i < uniformitySet.size(); i++){
			printf("[%s][%s]\n", uniformitySet[i][0].c_str(), uniformitySet[i][1].c_str());
		}
		printf("\n");
		printf("optiDifferSet 1213:\n");
		for (int i = 0; i < optiDifferSet.size(); i++){
			printf("[%s][%s]\n", optiDifferSet[i][0].c_str(), optiDifferSet[i][1].c_str());
		}
		printf("\n");

		printf("optiUniformitySet 1213:\n");
		for (int i = 0; i < optiUniformitySet.size(); i++){
			printf("[%s][%s]\n", optiUniformitySet[i][0].c_str(), optiUniformitySet[i][1].c_str());
		}
		printf("\n");
		*/
		
		//������֪����ȷ��˫��˫ȱ��
		//ȷ�����ԣ����ȸ���term���ɱ��壬�����LNF��صģ�������LNF��ص�ȱ��
		DnfMutant dnfMutant;
		//vector<vector<string>> termCharacters;//�������
		int shrinkNum = 0;//�����ĸ���
		int expandNum = 0;//���ŵĸ���
		int shrinkAndExpendNum = 0;//����+���ŵĸ���
		vector<string> maybeFaultTerms;

		for (int i = 0; i < termCharacters.size(); i++) {
			if (termCharacters[i][1] == "shrink"){
				shrinkNum++;
			}
			else if (termCharacters[i][1] == "expand"){
				expandNum++;
			}
			else{
				shrinkAndExpendNum++;
			}
			maybeFaultTerms.push_back(termCharacters[i][0]);
		}
	
		printf("shrinkNum: %d\n", shrinkNum);
		printf("expandNum: %d\n", expandNum);
		printf("shrinkAndExpendNum: %d\n", shrinkAndExpendNum);

		Mutants tempMutants(newExp);
		Mutants totalMutants(newExp);

		vector<string> maybeFaultTypes;
		//������
		if (termCharacters.size() == 2){

			//2�������� LIF/TOF   2x2
			if (shrinkNum == 2){
				maybeFaultTypes.push_back("LIFxLIF");	//1��LIFxLIF
				maybeFaultTypes.push_back("LIFxTOF");	//2��LIFxTOF
				maybeFaultTypes.push_back("TOFxTOF");	//3��TOFxTOF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2�����ţ� LOF/CORF  2x2
			else if (expandNum == 2){
				maybeFaultTypes.push_back("LOFxLOF");	//1��LOFxLOF
				maybeFaultTypes.push_back("LOFxCORF");	//2��LOFxCORF
				maybeFaultTypes.push_back("CORFxCORF");	//3��CORFxCORF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			// 1��������1�����ţ� LIF/TOF x LOF/CORF  2x2
			else if (shrinkNum == 1 && expandNum == 1){
				maybeFaultTypes.push_back("LIFxLOF");	//1��LIFxLOF
			
				maybeFaultTypes.push_back("TOFxLOF");	//3��TOFxLOF
				maybeFaultTypes.push_back("TOFxCORF");	//4��TOFxCORF

				maybeFaultTypes.push_back("LOFxLRF");	//1��LOFxLRF
				maybeFaultTypes.push_back("LOFxLNF");	//2��LOFxLNF
				maybeFaultTypes.push_back("LOFxTNF");	//3��LOFxTNF
				maybeFaultTypes.push_back("LIFxCORF");	//2��LIFxCORF
				maybeFaultTypes.push_back("CORFxLRF");	//4��CORFxLRF
				maybeFaultTypes.push_back("CORFxLNF");	//5��CORFxLNF
				maybeFaultTypes.push_back("CORFxTNF");	//6��CORFxTNF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//һ������ + һ������&���ţ� LIF/TOF + LRF/LNF/TNF  2 x 3
			else if (shrinkNum == 1 && shrinkAndExpendNum == 1){
				maybeFaultTypes.push_back("LIFxLRF");	//1��LIFxLRF
				maybeFaultTypes.push_back("LIFxLNF");	//2��LIFxLNF
				maybeFaultTypes.push_back("LIFxTNF");	//3��LIFxTNF
				maybeFaultTypes.push_back("TOFxLRF");	//4��TOFxLRF
				maybeFaultTypes.push_back("TOFxLNF");	//5��TOFxLNF
				maybeFaultTypes.push_back("TOFxTNF");	//6��TOFxTNF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//һ������ + һ������&���ţ� LOF/CORF + LRF/LNF/TNF  2 x 3
			else if (expandNum == 1 && shrinkAndExpendNum == 1){
				maybeFaultTypes.push_back("LOFxLRF");	//1��LOFxLRF
				maybeFaultTypes.push_back("LOFxLNF");	//2��LOFxLNF
				maybeFaultTypes.push_back("LOFxTNF");	//3��LOFxTNF
				maybeFaultTypes.push_back("CORFxLRF");	//4��CORFxLRF
				maybeFaultTypes.push_back("CORFxLNF");	//5��CORFxLNF
				maybeFaultTypes.push_back("CORFxTNF");	//6��CORFxTNF

				//maybeFaultTypes.push_back("CORFxLIF");	//����: CORFxLIF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//��������&���ţ� LRF/LNF/TNF + LRF/LNF/TNF  3 x 3
			else if (shrinkAndExpendNum == 2){
				maybeFaultTypes.push_back("LIFxLRF");	//2��LIFxLNF
				maybeFaultTypes.push_back("LRFxLRF");	//1��LRFxLRF
				maybeFaultTypes.push_back("LRFxLNF");	//2��LRFxLNF
				maybeFaultTypes.push_back("LRFxTNF");	//3��LRFxTNF
				maybeFaultTypes.push_back("LNFxLNF");	//4��LNFxLNF
				maybeFaultTypes.push_back("LNFxTNF");	//5��LNFxTNF
				maybeFaultTypes.push_back("TNFxTNF");	//6��TNFxTNF
			
				maybeFaultTypes.push_back("LIFxLNF");	//2��LIFxLNF
				maybeFaultTypes.push_back("LIFxTNF");	//3��LIFxTNF
				maybeFaultTypes.push_back("TOFxLRF");	//4��TOFxLRF
				maybeFaultTypes.push_back("TOFxLNF");	//5��TOFxLNF
				maybeFaultTypes.push_back("TOFxTNF");	//6��TOFxTNF


				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}
			
		}

		//������� �϶���DORF �� ENF
		else if (termCharacters.size() == 3 || termCharacters.size() >= 4){

			//2������ + 1������    DORF x LIF/TOF    1x2
			if (shrinkNum == 3){
				maybeFaultTypes.push_back("DORFxLIF");	//1��DORFxLIF
				maybeFaultTypes.push_back("DORFxTOF");	//2��DORFxTOF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2������ + 1������    DORF x LOF/CORF    1x2
			else if (shrinkNum == 2 && expandNum == 1){
				maybeFaultTypes.push_back("DORFxLOF");	//1��DORFxLOF
				maybeFaultTypes.push_back("DORFxCORF");	//2��DORFxCORF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2������&���� + 1������     ENF x LIF/TOF    1x2
			else if (shrinkAndExpendNum == 2 && shrinkNum == 1){
				maybeFaultTypes.push_back("ENFxLIF");	//1��ENFxLIF
				maybeFaultTypes.push_back("ENFxTOF");	//2��ENFxTOF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2������&���� + 1������   ENF x LOF/CORF    1x2
			else if (shrinkAndExpendNum == 2 && expandNum == 1){
				maybeFaultTypes.push_back("ENFxLOF");	//1��ENFxLOF
				maybeFaultTypes.push_back("ENFxCORF");	//2��ENFxCORF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}
		}
		maybeFaultTypes.push_back("DORFxTNF");	//1��DORFxTNF
		maybeFaultTypes.push_back("DORFxLRF");	//2��DORFxLRF
		maybeFaultTypes.push_back("DORFxDORF");	//3��DORFxDORF
		maybeFaultTypes.push_back("DORFxLNF");	//4��DORFxTNF
		maybeFaultTypes.push_back("DORFxLOF");	//4��DORFxLOF

		maybeFaultTypes.push_back("ENFxLRF");	//5��ENFxLRF
		maybeFaultTypes.push_back("ENFxLNF");	//6��ENFxLNF
		maybeFaultTypes.push_back("ENFxTNF");	//7��ENFxTNF
		maybeFaultTypes.push_back("ENFxDORF");	//8��ENFxDORF
		maybeFaultTypes.push_back("ENFxENF");	//9��ENFxENF


		maybeFaultTypes.push_back("TNFxTNF");	//10��TNFxTNF

		maybeFaultTypes.push_back("TOFxTNF");	//11��TOFxTNF

		maybeFaultTypes.push_back("TOFxCORF");	//12��TOFxCORF
		maybeFaultTypes.push_back("CORFxCORF");	//13��CORFxCORF

		maybeFaultTypes.push_back("CORFxLNF");	//14��CORFxLNF
		maybeFaultTypes.push_back("CORFxLOF");	//14��CORFxLOF
		maybeFaultTypes.push_back("CORFxLIF");	//14��CORFxLIF
		maybeFaultTypes.push_back("CORFxLRF");	//14��CORFxLRF

		for (int i = 0; i < maybeFaultTypes.size(); i++){
			string faultType = maybeFaultTypes.at(i);
			if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
				return 0;
			}
		}
	}
	return 0;
}

