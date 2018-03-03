// BoolCheckFramework.cpp : 定义控制台应用程序的入口点。
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

//根据表达式获取变体
//参数exp，为原表达式
//参数faultType，为缺陷类型
//参数iModel，为范式类型
//生成DNF表达式
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
	vector<string> testset/*生成的测试用例集合，如果测试用例全部通过，则表示没有问题*/,
	string exprOri/*原表达式*/,
	string expr/*变体*/,
	vector<int>& unformityTCsol, /*结果一致的节点，为无效测试点*/
	vector<int>& differTCsol, /*结果不一致的节点，为有效测试点*/
	vector<string>& unformityStr, /*结果一致的节点，为无效测试点*/
	vector<string>& differTStr, /*结果不一致的节点，为有效测试点*/
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
		bool v = v1 ^ v2;//两个值相同则为假，不同则为真,不过这两个资源好浪费时间呀，明明有更优的算法
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

//剔除约束矩阵中全零的行
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
		printf("表达式只包含一项，不可能发生这种Double Fault LIFxLIF！");
		return false;
	}
	//double fault lifxlif 每次需要处理两项
	// 1.1预处理，将每个term的dimensionValuesStr存起来
	vector<string> dimensionValuesStrs; //存储每个term的dimension value string
	int termNum = terms.size();
	for (int i = 0; i<termNum; i++)//逐个term进行解析
	{
		string dimensionValuesStr = "";
		for (int j = 0; j < dimensionNum; j++)
		{
			//将各维的01表示法，转换为1,2,3,4表示法
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			if (j != dimensionNum - 1)
				dimensionValuesStr.append("\n");
		}
		dimensionValuesStrs.push_back(dimensionValuesStr);
	}
	//1.2 预处理读取三项以上otp
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
		printf("表达式只包含一项，不可能发生这种Double Fault LIFxLRF！");
		return false;
	}
	//double fault lifxlrf 每次需要处理两项
	// 1.1预处理，将每个term的dimensionValuesStr存起来
	int termNum = terms.size();
	vector<string> dimensionValuesStrs; //存储每个term的dimension value string
	for (int i = 0; i<termNum; i++)//逐个term进行解析
	{
		string dimensionValuesStr = "";
		for (int j = 0; j < dimensionNum; j++)
		{
			//将各维的01表示法，转换为1,2,3,4表示法
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			if (j != dimensionNum - 1)
				dimensionValuesStr.append("\n");
		}
		dimensionValuesStrs.push_back(dimensionValuesStr);
	}
	//1.2 预处理读取三项以上otp
	vector<int> all3OPs;
	for (int i = 0; i<opover2.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(opover2.at(i), vOdd);
		all3OPs.push_back(value);
	}
	sort(all3OPs.begin(), all3OPs.end());

	//生成LRF矩阵
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
		printf("表达式只包含一项，不可能发生这种Double Fault LRFxLRF！");
		return false;
	}
	//double fault lrfxlrf 每次需要处理两项
	// 1.1预处理，将每个term的dimensionValuesStr存起来
	vector<string> dimensionValuesStrs; //存储每个term的dimension value string
	for (int i = 0; i<termNum; i++)//逐个term进行解析
	{
		string dimensionValuesStr = "";
		bool isSingleLitral = false;
		//格子占全的维度数目
		int fullDemNum = 0;
		//格子占半的维度数目
		int halfDemNum = 0;
		for (int j = 0; j < dimensionNum; j++)
		{
			//将各维的01表示法，转换为1,2,3,4表示法
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			int strLength = tempStr.size();
			if (vOdd != true || j != dimensionNum - 1) //格子占全时候字符串长度应为7,占半长度为3
			{
				if (strLength == 7)
					fullDemNum++;
				else if (strLength == 3)
					halfDemNum++;
			}
			else //格子占全时候字符串长度应为3,占半长度为1
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
		//生成DLRF矩阵
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
			if (!v1)continue;//不是收缩点，舍去
			bool v2 = uti.evaluateMutant(oriExp, terms[i], testcase[j]);
			bool v = v1 ^ v2;//异或，值同为假，值不同为真
			if (v == false){//如果为假，说明这个项发生了缺陷
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
	//分解表达式
	PRE_PROCESS_DATA preData;

	//转换为IDNF
	Utility uti;
	vector<string> terms = uti.extractDNFstyle(expression);
	vector<string>::iterator iter;
	string newExp = "";
	//得到IDNF
	newExp = expression;// generateIDNForICNF(expression, true);
	//printf("newExp:%s\n", newExp.c_str());

	getPreProcessData(newExp, preData, DNF);
	CCommonFuncs commonFuncs;
	int literalNum = model.size();
	int dimensionNum = (literalNum + 1) / 2;
	bool vOdd = (literalNum % 2 == 0) ? false : true;
	//printf("%d,%d,%d\n", __LINE__, vOdd, literalNum);

	//获取变化词的维度编号
	int literalNo1;
	int dimensionNo1 = getDimensionNo(model, Liter1, literalNo1);

	//读取所有的OTP
	vector<string> op = preData.op;
	vector<int> allOPs;
	for (int i = 0; i<op.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(op.at(i), vOdd);
		allOPs.push_back(value);
	}
	sort(allOPs.begin(), allOPs.end());

	//读取所有的UTP
	vector<int> allUPs;
	vector<string> up = preData.up;
	for (int i = 0; i<up.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(up.at(i), vOdd);
		allUPs.push_back(value);
	}
	sort(allUPs.begin(), allUPs.end());

	//获取是否有单词的判断
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
	//一个term对应的维度值，有几个维度，dimensionNum就会有多大
	vector<string> testcase;
	CCommonFuncs commonFunc;
	mutantChangeTerm = mutantTerm;
	string lm = Liter1;
	int lz = mutantChangeTerm.find(lm);
	string::iterator itr = mutantChangeTerm.begin() + lz;
	mutantChangeTerm.erase(itr);
	//printf("LOF:mutantChangeTerm[%s]\n", mutantChangeTerm.c_str());

	vector<string> addValues;
	//首先判断是哪个维度，其次判断是
	/*printf("termsDimension[%d]:\n", dimensionNo1);
	for (int x = 0; x < termsDimension[dimensionNo1].size(); x++){
		printf("[%s]", termsDimension[dimensionNo1][x].c_str());
	}
	printf("\n");*/
	//没有对非做处理
	if (literalNo1 == dimensionNo1 * 2){//如果是4个数，则不变，否则增加该词的01值
		//新增，如果是1，则增加为2，如果是2，则增加为4
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
	//变化后的格覆盖
	/*printf("\nchanged:\n");
	for (int x = 0; x < dimensionNum; x++){
		for (int y = 0; y < termsDimension[x].size(); y++){
			printf("[%s]", termsDimension[x][y].c_str());
		}
		printf("\n");
	}*/
	//变化的格覆盖点
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


bool checkTestDiffer(string test, string exprOri/*原表达式*/, string expr/*变体*/){
	Utility uti;
	bool v = (uti.evaluate(exprOri, test)) ^ (uti.evaluateMutant(exprOri, expr, test));//异或：不同为真，相同为假
	return v;
}

//诊断——用最优测试用例
bool diagnosis(
	string oriExp, //原表达式
	string faultExp, //待测表达式
	string faultType, //缺陷类型
	vector<string> oriTerms,
	vector<string> faultTerms,
	string model,
	vector<vector<string>> optiUniformitySet,//最优无效点
	vector<vector<string>> optiDifferSet//最优有效点
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
		bool flagUniform = true;//无效点是否全部通过
		//（1）先用无效点进行探测
		for (int i = 0; i < optiUniformitySet.size(); i++){
			//printf("point: %s, faultExp: %s, mutant: %s\n", optiUniformitySet[i][0].c_str(), faultExp.c_str(), mutant.c_str());
			if (checkTestDiffer(optiUniformitySet[i][0], faultExp, mutant)){
				flagUniform = false;
				break;
			}
		}

		if (flagUniform == true){//无效点全部通过
			bool flagDiffer = true;
			for (int i = 0; i < optiDifferSet.size(); i++){
				if (checkTestDiffer(optiDifferSet[i][0], faultExp, mutant)){
					flagDiffer = false;
					break;
				}
			}
			if (flagDiffer == true){//有效点全部通过
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
	//分解表达式
	PRE_PROCESS_DATA preData;

	clock_t time1 = clock();
	//转换为IDNF
	Utility uti;
	vector<string> terms = uti.extractDNFstyle(expression);
	string model = uti.generateModel(expression);
	string newExp = "";
	//得到IDNF
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
	//读取所有的OTP
	vector<string> op = preData.op;
	vector<int> allOPs;
	for (int i = 0; i<op.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(op.at(i), vOdd);
		allOPs.push_back(value);
	}
	sort(allOPs.begin(), allOPs.end());

	//读取所有的UTP·
	vector<int> allUPs;
	vector<string> up = preData.up;
	for (int i = 0; i<up.size(); i++)
	{
		int value = commonFuncs.getTModelIndexOf01Str(up.at(i), vOdd);
		allUPs.push_back(value);
	}
	sort(allUPs.begin(), allUPs.end());

	//获取是否有单词的判断
	bool isSingleLitral = false;
	//格子占全的维度数目
	int fullDemNum = 0;
	//格子占半的维度数目
	int halfDemNum = 0;
	vector<vector<string>> termsDimension = preData.termsDimension;

	vector<vector<vector<int>>> mutantsDimensionValue;
	vector<vector<CCell*>> corresponding_tests;
	vector<int> tModelIndexAppeared;
	vector<vector<CCell*>> matrixCells;
	vector<vector<int>> coefficientMatrix;
	vector<string> mutants;
	/*printf("%s,%d,termsDimension Print:\n", __FILE__, __LINE__);//得到该term的维度变化，一行为一个term的变化维度
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
	for (int i = 0; i < preData.termsDimension.size(); i++){//逐个term进行解析
		string dimensionValuesStr;
		//将一行转换为多维变换，例如11 10 00 01 11 10 ，就是3 4 1 2 3 4，即34 1234
		for (int j = 0; j < dimensionNum; j++)
		{
			//将各维的01表示法，转换为1,2,3,4表示法
			string line = termsDimension[i][j];
			string tempStr;
			commonFuncs.transDimensionLineFrom01totModelIndex(line, tempStr);
			dimensionValuesStr.append(tempStr);
			int strLength = tempStr.length();
			if (vOdd != true || j != dimensionNum - 1) //格子占全时候字符串长度应为7,占半长度为3
			{
				if (strLength == 7)
					fullDemNum++;
				else if (strLength == 3)
					halfDemNum++;
			}
			else //格子占全时候字符串长度应为3,占半长度为1
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
		//根据缺陷信息转换成需要的格覆盖,在此处要进行改造
		mutantsDimensionValue.clear();
		corresponding_tests.clear();
		//整个思路
		//首先根据本体得到所有变体的待测点，所谓待测点，就是变体与本体不同的真值点，即包括本体真值点，也包括变体真值点
		//然后根据待测点组方程组，一个变体一个方程组，方程组变量的个数为所有待测点
		//然后根据每个变体的待测点确定每个方程的变量个数
		//结方程组得到最优的变量个数即可，其中不涉及到表达式的真假值，至于整个处理的逻辑原理请参考表达式
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
			if (isSingleLitral) //如果只有一个词 则退化成TOF
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

			//添加6个非CC-S可检测的双缺陷
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

	//转换约束条件
	//得到所有的真值点，这一步必须进入循环吗？
	commonFuncs.getAlltModelIndexAppearedInMatrix(matrixCells, tModelIndexAppeared);
	//对列向量进行排序，如311排在312之前
	sort(tModelIndexAppeared.begin(), tModelIndexAppeared.end());

	commonFuncs.getCoefficientMatrix(tModelIndexAppeared, matrixCells, coefficientMatrix);

	vector<vector<int>>finalCoefficientMatrix;
	vector<int>finalTModelIndex;
	removeAllZeroCols(coefficientMatrix, tModelIndexAppeared, finalCoefficientMatrix, finalTModelIndex);

	vector<int>().swap(tModelIndexAppeared);//只是为了打印用
	vector<vector<int>>().swap(coefficientMatrix);//只是为了打印用
	matrixCells.clear();
	//Glpk解约束方程
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
	//可以修改，不需要每次都生成
	vector<string> mutants_DP;
	for (int num = 0; num<mutants.size(); num++)
	{
		//printf("%s,%d,mutants[%d]:[%s]\r\n", __FILE__, __LINE__, num, mutants[num].c_str());
		mutants_DP.push_back(mutants[num]);
	}

	glpk_Solve.run(colSelected, &finalCoefficientMatrix);
	//将得到的格覆盖转换成测试用例
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
	//printf("测试用例生成完毕[%s]\n", faultExpression.c_str());

	vector<vector<int>> DP_Matrix;
	vector<vector<int>> case_mutant;

	genDP_Matrix(finalCoefficientMatrix, finalTModelIndex, DP_Matrix);//对约束矩阵进行置转，这样一行就对应一个格覆盖

	genCase_Mutants(case_mutant, finalTModelIndex, DP_Matrix);
	vector<int> finalTModelIndexMatrixCPY(finalTModelIndex);
	vector<vector<int>> finalCoefficientMatrixCPY(finalCoefficientMatrix);
	FaultDetection faultDetector;
	if (faultExpression == ""){
		//检测表达式
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
		{//不用循环找，直接找方程即可
			vector<int> uniformityTCsol, differTCsol;
			vector<string> uniformityStr, differStr;
			vector<string> uniformityType, differType;
			//获取有效点和无效点
			if (checkTestSet(
				testset,
				newExp, faultExpression,
				uniformityTCsol, differTCsol, uniformityStr, differStr, uniformityType, differType, vOdd))
			{//返回真，说明全是无效测试点
				printf("\ncannot find any faults\n");
				return 0;
			}
			else{
				bool isExpand = false;
				bool isShrink = false;
				//记录无效点
				for (int l = 0; l < uniformityStr.size(); l++){
					vector<string> llm;
					llm.push_back(uniformityStr[l]);
					llm.push_back(uniformityType[l]);
					uniformitySet.push_back(llm);
					if (isFirstO){
						optiUniformitySet.push_back(llm);
					}
				}
				//记录有效点
				for (int l = 0; l < differStr.size(); l++){
					if (differType[l] == "e"){//不满足条件
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

				if (uniformityTCsol.size() != 0)sort(uniformityTCsol.begin(), uniformityTCsol.end());//一致的点
				if (differTCsol.size() != 0)sort(differTCsol.begin(), differTCsol.end());//不一致的点
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
				glpk_Solve.run_lo(colSelected1, &DP_Matrix);//这里其实可以删除其他列再置转求解，结果是一样的，但老师没说让这样，没有理论就没有伤害
				/*printf("colSelected1 1211:\n");
				for (int i = 0; i < colSelected1.size(); i++){
					printf("[%d]", colSelected1[i]);
				}
				printf("\n");*/
				for (int j = 0; j < colSelected1.size(); j++)
				{
					//printf("mutant[%d][%s][%s]\n", colSelected1[j], mutants[colSelected1[j]].c_str(), faultTypes[colSelected1[j]].c_str());
					//如果是LIF，则特殊处理，只有收缩点则，该表达式保留，并保留对应的LRF，如果有扩张点，则不保留该表达式，并看对应的LRF，只有扩张点，则直接看LOF
					if (faultTypes[colSelected1[j]] == "LIF"){
						bool remainSign = false;
						bool tmt = true;
						if (isExpand){//只有扩张就不比较了，直接下一步
							tmt = false;
						}
						else{
							for (int y = 0; y < uniformityStr.size(); y++){//无效点满足条件
								bool v2 = (uti.evaluateMutant(newExp, mutants[colSelected1[j]], uniformityStr[y]));
								if (uti.evaluate(newExp, uniformityStr[y]) ^ v2){//不满足条件
									tmt = false;
									break;
								}
							}
						}
						if (tmt){
							//printf("LIF remain\n");
							remainSign = true;//所有点都正确，保留
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
						{//找LRF
							//直接找对应的LRF，生成方法，直接LOF，并剔除等价表达式
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
							Mutants ab = dnfMutant.generateLOFmutants(mutants[colSelected1[j]], faultOccurTerms);//LIF&LOF就是LRF
							for (int y = 0; y < ab.mutants.size(); y++){
								//如果该表达式与原表达式等价，则剔除
								if (checkEqualment(newExp, ab.mutants[y].mutant)){
									continue;
								}
								//printf("Find LRF mutant[%s][%s][%s]\n", ab.mutants[y].mutant.c_str(), "LRF", faultTerms[colSelected1[j]].c_str());
								//将测试用例带入检验是否正确
								bool tmt = true;
								for (int z = 0; z < uniformityStr.size(); z++){
									bool v2 = (uti.evaluateMutant(newExp, ab.mutants[y].mutant, uniformityStr[z]));
									if (uti.evaluate(newExp, uniformityStr[z]) ^ v2){//不满足条件
										tmt = false;
										break;
									}
								}
								if (tmt){//所有点都正确，保留
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

						if (mutants[colSelected1[j]] == ""){//没有生成变体，则直接检验约束方程
							/*不再考虑了，直接用单缺陷推测双缺陷*/
							int uvl = finalCoefficientMatrix[colSelected1[j]].size();
							int uvn = finalTModelIndex.size();
							bool rightSign = true;
							for (int x = 0; x < uniformityTCsol.size(); x++){
								for (int y = 0; y < finalTModelIndex.size(); y++){
									if (uniformityTCsol[x] == finalTModelIndex[y]){
										if (finalCoefficientMatrix[colSelected1[j]][y] != 0){//不满足条件
											rightSign = false;
											break;
										}
									}
								}
								if (!rightSign)break;
							}
							if (!rightSign)continue;
							else{//找到了一个约束方程
								//printf("maybe a mutant[%d][%s][%s][%s]\n", colSelected1[j], mutants[colSelected1[j]].c_str(), faultTypes[colSelected1[j]].c_str(), faultTerms[colSelected1[j]].c_str());
								if (find(selectAddRows.begin(), selectAddRows.end(), colSelected1[j]) == selectAddRows.end()){
									vector<string> ltsd;
									char lz[10];
									sprintf(lz, "#%d", colSelected1[j]);
									ltsd.push_back(lz);//直接记录标号好了
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
							{//返回真，说明全是无效测试点
								continue;
							}
							else{
								bool rightSign = true;
								if (faultTypes[colSelected1[j]] == "LOF"){
									//如果有收缩点，直接排除
									if (isShrink){
										rightSign = false;
									}
									else{
										for (int x = 0; x < uniformityTCsol.size(); x++){
											bool v2 = (uti.evaluateMutant(newExp, mutants[colSelected1[j]], uniformityStr[x]));
											if (uti.evaluate(newExp, uniformityStr[x]) ^ v2){//不满足条件
												rightSign = false;
												break;
											}
										}
									}
								}
								else{
									for (int x = 0; x < uniformityTCsol.size(); x++){
										bool v2 = (uti.evaluateMutant(newExp, mutants[colSelected1[j]], uniformityStr[x]));
										if (uti.evaluate(newExp, uniformityStr[x]) ^ v2){//不满足条件
											rightSign = false;
											break;
										}
									}
								}
								if (!rightSign)continue;
								else{//找到了一个约束方程
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
				//将所有满足条件的LOF和其他类型的变体或方程都记录在案，使用无效点判断
				for (int x = 0; x < finalCoefficientMatrix.size(); x++){
					if (find(selectAddRows.begin(), selectAddRows.end(), x) == selectAddRows.end()){//说明没有该记录
						bool ntl = true;
						for (int i = 0; i < uniformityTCsol.size() && ntl; i++){
							for (int j = 0; j < finalTModelIndex.size() && ntl; j++){
								if (uniformityTCsol[i] == finalTModelIndex[j]){//找到了这个变量，看看哪些方程可以挂钩
									if (finalCoefficientMatrix[x][j] != 0){//无效点
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
		//检查是否还有没有验证的测试用例，如果有，则直接进行检验，用于验证变体
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
		//获取有效点和无效点
		checkTestSet(
			testset,
			newExp, faultExpression,
			uniformityTCsol, differTCsol, uniformityStr, differStr, uniformityType, differType, vOdd);

		//这里测试根据有收缩点判断扩张点对应的term
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
		//首先确定是否有扩张点，如果有，确定该项的所有扩张点
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

		for (int i = shrinkTermsNo.size() - 1; i >= 0; i--){//判断收缩点中是否有扩张点
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

		//记录无效点
		for (int l = 0; l < uniformityStr.size(); l++){
			vector<string> llm;
			llm.push_back(uniformityStr[l]);
			llm.push_back(uniformityType[l]);
			uniformitySet.push_back(llm);
		}
		//记录有效点
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

		//使用有效点的收缩扩张过滤maybeMutants
		bool sSign = false;
		bool eSign = false;
		for (int j = 0; j < differSet.size(); j++){
			if (differSet[j][1] == "s")sSign = true;
			if (differSet[j][1] == "e")eSign = true;
			if (sSign && eSign)break;
		}
		for (int i = maybeMutants.size() - 1; i >= 0;){
			bool lmt = false;
			if (sSign && !eSign){//只有收缩，剔除LOF
				if (maybeMutants[i][1] == "LOF" || maybeMutants[i][1] == "LRF"){//此时LRF有LIF的等价表达式
					vector<vector<string>>::iterator itr = maybeMutants.begin() + i;
					maybeMutants.erase(itr);
					i--;
					if (i < 0)break;
					lmt = true;
				}
			}
			if (eSign && !sSign){//只有扩张，剔除LIF
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
		//使用无效点正确性过滤maybeMutants
		for (int i = maybeMutants.size() - 1; i >= 0; i--){
			bool selectedSign = true;//新逻辑中如果是双缺陷，不再考虑没有变体的方程
			if (maybeMutants[i][0][0] == '#'){//没有变体
				string l = maybeMutants[i][0].substr(1, maybeMutants[i][0].size() - 1);
				int ln = atoi(l.c_str());
				for (int i = 0; i < uniformityTCsol.size() && selectedSign; i++){
					for (int j = 0; j < finalTModelIndex.size() && selectedSign; j++){
						if (uniformityTCsol[i] == finalTModelIndex[j]){//找到了这个变量，看看哪些方程可以挂钩
							if (finalCoefficientMatrix[ln][j] != 0){//无效点
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
					if (v1 ^ v2){//只要有一个不对就剔除
						selectedSign = false;
						break;
					}
				}
				//使用变体的无效点过滤
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
				if (!(v1 ^ v2)){//只要有一个不对就剔除
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

		//如果maybeMutants中没有数据，说明肯定是双缺陷
		//如果有数据且为LIF或LOF有两个term，表示肯定是不含LRF的双向双缺陷
		//如果有数据且为LIF或LOF有一个term，表示肯定为单缺陷，或不含LRF的单项双缺陷
		//如果有数据且为LRF有一个term，表示肯定为单缺陷，或含LRF的单项双缺陷
		//如果没有数据，则采用反项无效点过滤,即有效点过滤
		if (maybeMutants.size() == 0){
			for (int i = maybeMutants.size() - 1; i >= 0; i--){
				bool selectedSign = true;
				for (int j = 0; j < differSet.size(); j++){
					bool v1 = (uti.evaluate(newExp, differSet[j][0]));
					bool v2 = (uti.evaluateMutant(newExp, maybeMutants[i][0], differSet[j][0]));
					if (v1 ^ v2){//只要有一个不对就剔除
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

		//如果表达式是我们需要的表达式，则直接退出：
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
		//将缺陷类型和发生缺陷项进行过滤
		vector<string>maybeTypes;
		vector<vector<string>>maybeTerms;
		//归并缺陷类型
		for (int i = 0; i < maybeMutants.size(); i++){
			bool lmdt = true;
			for (int j = 0; j < maybeTypes.size(); j++){
				if (maybeMutants[i][1] == maybeTypes[j]){//找到了缺陷类型则不用插入
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
		//归并缺陷发生项
		for (int i = 0; i < maybeMutants.size(); i++){
			bool lmdt = true;
			for (int j = 0; j < maybeTypes.size(); j++){
				if (maybeMutants[i][1] == maybeTypes[j]){//找到了缺陷类型
					bool muv = true;
					for (int k = 0; k < maybeTerms[j].size(); k++){
						if (maybeTerms[j][k] == maybeMutants[i][2]){//该项已经存在
							muv = false;
						}
					}
					if (muv){//该项不存在
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

		//如果类型中没有LIF，直接剔除LRF，如果只有一个LIF，直接剔除LRFxLRF和有多个项中的LRF中的其他项
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
			//剔除LRF和LRFxLRF
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
				//剔除LRFxLRF和LRF多项中的其他项
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

		//后续接口打印：
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
		
		//根据已知条件确定双向双缺陷
		//确定策略，首先根据term生成变体，如果是LNF相关的，则生成LNF相关的缺陷
		DnfMutant dnfMutant;
		//vector<vector<string>> termCharacters;//项的特征
		int shrinkNum = 0;//收缩的个数
		int expandNum = 0;//扩张的个数
		int shrinkAndExpendNum = 0;//收缩+扩张的个数
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
		//两个项
		if (termCharacters.size() == 2){

			//2个收缩， LIF/TOF   2x2
			if (shrinkNum == 2){
				maybeFaultTypes.push_back("LIFxLIF");	//1、LIFxLIF
				maybeFaultTypes.push_back("LIFxTOF");	//2、LIFxTOF
				maybeFaultTypes.push_back("TOFxTOF");	//3、TOFxTOF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2个扩张， LOF/CORF  2x2
			else if (expandNum == 2){
				maybeFaultTypes.push_back("LOFxLOF");	//1、LOFxLOF
				maybeFaultTypes.push_back("LOFxCORF");	//2、LOFxCORF
				maybeFaultTypes.push_back("CORFxCORF");	//3、CORFxCORF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			// 1个收缩，1个扩张， LIF/TOF x LOF/CORF  2x2
			else if (shrinkNum == 1 && expandNum == 1){
				maybeFaultTypes.push_back("LIFxLOF");	//1、LIFxLOF
			
				maybeFaultTypes.push_back("TOFxLOF");	//3、TOFxLOF
				maybeFaultTypes.push_back("TOFxCORF");	//4、TOFxCORF

				maybeFaultTypes.push_back("LOFxLRF");	//1、LOFxLRF
				maybeFaultTypes.push_back("LOFxLNF");	//2、LOFxLNF
				maybeFaultTypes.push_back("LOFxTNF");	//3、LOFxTNF
				maybeFaultTypes.push_back("LIFxCORF");	//2、LIFxCORF
				maybeFaultTypes.push_back("CORFxLRF");	//4、CORFxLRF
				maybeFaultTypes.push_back("CORFxLNF");	//5、CORFxLNF
				maybeFaultTypes.push_back("CORFxTNF");	//6、CORFxTNF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//一个收缩 + 一个收缩&扩张， LIF/TOF + LRF/LNF/TNF  2 x 3
			else if (shrinkNum == 1 && shrinkAndExpendNum == 1){
				maybeFaultTypes.push_back("LIFxLRF");	//1、LIFxLRF
				maybeFaultTypes.push_back("LIFxLNF");	//2、LIFxLNF
				maybeFaultTypes.push_back("LIFxTNF");	//3、LIFxTNF
				maybeFaultTypes.push_back("TOFxLRF");	//4、TOFxLRF
				maybeFaultTypes.push_back("TOFxLNF");	//5、TOFxLNF
				maybeFaultTypes.push_back("TOFxTNF");	//6、TOFxTNF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//一个扩张 + 一个收缩&扩张， LOF/CORF + LRF/LNF/TNF  2 x 3
			else if (expandNum == 1 && shrinkAndExpendNum == 1){
				maybeFaultTypes.push_back("LOFxLRF");	//1、LOFxLRF
				maybeFaultTypes.push_back("LOFxLNF");	//2、LOFxLNF
				maybeFaultTypes.push_back("LOFxTNF");	//3、LOFxTNF
				maybeFaultTypes.push_back("CORFxLRF");	//4、CORFxLRF
				maybeFaultTypes.push_back("CORFxLNF");	//5、CORFxLNF
				maybeFaultTypes.push_back("CORFxTNF");	//6、CORFxTNF

				//maybeFaultTypes.push_back("CORFxLIF");	//特殊: CORFxLIF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//两个收缩&扩张， LRF/LNF/TNF + LRF/LNF/TNF  3 x 3
			else if (shrinkAndExpendNum == 2){
				maybeFaultTypes.push_back("LIFxLRF");	//2、LIFxLNF
				maybeFaultTypes.push_back("LRFxLRF");	//1、LRFxLRF
				maybeFaultTypes.push_back("LRFxLNF");	//2、LRFxLNF
				maybeFaultTypes.push_back("LRFxTNF");	//3、LRFxTNF
				maybeFaultTypes.push_back("LNFxLNF");	//4、LNFxLNF
				maybeFaultTypes.push_back("LNFxTNF");	//5、LNFxTNF
				maybeFaultTypes.push_back("TNFxTNF");	//6、TNFxTNF
			
				maybeFaultTypes.push_back("LIFxLNF");	//2、LIFxLNF
				maybeFaultTypes.push_back("LIFxTNF");	//3、LIFxTNF
				maybeFaultTypes.push_back("TOFxLRF");	//4、TOFxLRF
				maybeFaultTypes.push_back("TOFxLNF");	//5、TOFxLNF
				maybeFaultTypes.push_back("TOFxTNF");	//6、TOFxTNF


				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}
			
		}

		//三个项—— 肯定有DORF 或 ENF
		else if (termCharacters.size() == 3 || termCharacters.size() >= 4){

			//2个收缩 + 1个收缩    DORF x LIF/TOF    1x2
			if (shrinkNum == 3){
				maybeFaultTypes.push_back("DORFxLIF");	//1、DORFxLIF
				maybeFaultTypes.push_back("DORFxTOF");	//2、DORFxTOF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2个收缩 + 1个扩张    DORF x LOF/CORF    1x2
			else if (shrinkNum == 2 && expandNum == 1){
				maybeFaultTypes.push_back("DORFxLOF");	//1、DORFxLOF
				maybeFaultTypes.push_back("DORFxCORF");	//2、DORFxCORF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2个收缩&扩张 + 1个扩张     ENF x LIF/TOF    1x2
			else if (shrinkAndExpendNum == 2 && shrinkNum == 1){
				maybeFaultTypes.push_back("ENFxLIF");	//1、ENFxLIF
				maybeFaultTypes.push_back("ENFxTOF");	//2、ENFxTOF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}

			//2个收缩&扩张 + 1个扩张   ENF x LOF/CORF    1x2
			else if (shrinkAndExpendNum == 2 && expandNum == 1){
				maybeFaultTypes.push_back("ENFxLOF");	//1、ENFxLOF
				maybeFaultTypes.push_back("ENFxCORF");	//2、ENFxCORF

				for (int i = 0; i < maybeFaultTypes.size(); i++){
					string faultType = maybeFaultTypes.at(i);
					if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
						return 0;
					}
				}
			}
		}
		maybeFaultTypes.push_back("DORFxTNF");	//1、DORFxTNF
		maybeFaultTypes.push_back("DORFxLRF");	//2、DORFxLRF
		maybeFaultTypes.push_back("DORFxDORF");	//3、DORFxDORF
		maybeFaultTypes.push_back("DORFxLNF");	//4、DORFxTNF
		maybeFaultTypes.push_back("DORFxLOF");	//4、DORFxLOF

		maybeFaultTypes.push_back("ENFxLRF");	//5、ENFxLRF
		maybeFaultTypes.push_back("ENFxLNF");	//6、ENFxLNF
		maybeFaultTypes.push_back("ENFxTNF");	//7、ENFxTNF
		maybeFaultTypes.push_back("ENFxDORF");	//8、ENFxDORF
		maybeFaultTypes.push_back("ENFxENF");	//9、ENFxENF


		maybeFaultTypes.push_back("TNFxTNF");	//10、TNFxTNF

		maybeFaultTypes.push_back("TOFxTNF");	//11、TOFxTNF

		maybeFaultTypes.push_back("TOFxCORF");	//12、TOFxCORF
		maybeFaultTypes.push_back("CORFxCORF");	//13、CORFxCORF

		maybeFaultTypes.push_back("CORFxLNF");	//14、CORFxLNF
		maybeFaultTypes.push_back("CORFxLOF");	//14、CORFxLOF
		maybeFaultTypes.push_back("CORFxLIF");	//14、CORFxLIF
		maybeFaultTypes.push_back("CORFxLRF");	//14、CORFxLRF

		for (int i = 0; i < maybeFaultTypes.size(); i++){
			string faultType = maybeFaultTypes.at(i);
			if (diagnosis(newExp, faultExpression, faultType, terms, maybeFaultTerms, model, uniformitySet, differSet)){
				return 0;
			}
		}
	}
	return 0;
}

