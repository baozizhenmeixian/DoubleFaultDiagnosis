#include "stdafx.h"
#include "Utility.h"
#include<algorithm>

Utility::Utility(void)
{
}


Utility::~Utility(void)
{
}

string Utility::getLiteralNegation(string literal)
{
	if (literal.at(0) == '!')
	{
		literal.erase(0, 1);//ɾ��������
	}
	else
	{
		literal = "!" + literal;
	}
	return literal;
}

//��term�����ڲ����ʵ�����b!ac --> !abc
string Utility::sortTerm(string term)
{
	string newTerm;
	vector<string> literalTerm;
	vector<int> orderIndex;

	for (int i = 0; i < term.length(); i++)
	{
		char ch = term[i];
		if (ch != '!'){
			int num = ch - 'a';
			orderIndex.push_back(num);
			sort(orderIndex.begin(), orderIndex.end());
			int index = 0;
			for (int j = 0; j < orderIndex.size(); j++){
				if (orderIndex[j] == num){
					index = j;
					break;
				}
			}
			stringstream stream;
			stream << ch;
			string str = stream.str();
			literalTerm.insert(literalTerm.begin() + index, str);
		}
		else{
			ch = term[i + 1];
			i++;
			int num = ch - 'a';
			orderIndex.push_back(num);
			sort(orderIndex.begin(), orderIndex.end());
			int index = 0;
			for (int j = 0; j < orderIndex.size(); j++){
				if (orderIndex[j] == num){
					index = j;
					break;
				}
			}
			stringstream stream;
			stream << ch;
			string str = "!" + stream.str();
			literalTerm.insert(literalTerm.begin() + index, str);
		}
	}
	for (int i = 0; i < literalTerm.size(); i++){
		newTerm += literalTerm[i];
	}
	return newTerm;
}

vector<string> Utility::uniqueTermVector(vector<string> termVector, string term)
{
	bool flag = true;
	for (int i = 0; i < termVector.size(); i++){
		if (termVector[i] == term){
			flag = false;
		}
	}
	if (flag){
		termVector.push_back(term);
	}
	return termVector;
}
//For example, extract "a!d(e(b+c)+a)" into three terms: "a" ,"!d" and "(e(b+c)+a)"
vector<string> Utility::extractCNFstyle(string exp)
{
	vector<string> result;

	int inBrackets=0;
	int bracketTermStart=0;
	for(int i=0;i<exp.length();i++)
	{
		wchar_t ch=exp[i];
		if(inBrackets==0){
			if(ch=='('){
				inBrackets++;
				bracketTermStart=i;
			}
			else if(ch=='!')
				continue;
			else{
				string s;
				if(i!=0&&exp[i-1]=='!'){
					s.push_back('!');
					s.push_back(ch);
					result.push_back(s);
				}
				else{
					s.push_back(ch);
					result.push_back(s);
				}
			}
		}else{
			if(ch==')'&&inBrackets==1){
				if(bracketTermStart!=0&&exp[bracketTermStart-1]=='!')
					result.push_back(exp.substr(bracketTermStart-1,i+1-(bracketTermStart-1)));
				else
					result.push_back(exp.substr(bracketTermStart,i+1-bracketTermStart));
				inBrackets--;
			}else if(ch=='(')
				inBrackets++;
			else if(ch==')')
				inBrackets--;
			else 
				continue;
		}
	}//for
	return result;
}

//For example, extract "a!d(e(b+c)+a)+c(d+e)" into two terms: "a!d(e(b+c)+a)" and "c(d+e)"
vector<string> Utility::extractDNFstyle(string exp)
{
	vector<string> dnfStyle;
	int inBrackets=0;
	int termStart=0;
	for(int i=0;i<exp.length();i++)
	{
		wchar_t ch=exp[i];
		if(ch=='(')
			inBrackets++;
		else if(ch==')')
			inBrackets--;
		else if(ch=='+'&&inBrackets==0){
			string term=exp.substr(termStart,i-termStart);
			dnfStyle.push_back(term);
			termStart=i+1;
		}
	}//for
	string term=exp.substr(termStart);
	dnfStyle.push_back(term);
	return dnfStyle;
}

//result will be sorted as the sequence of alphabet
string Utility::generateModel(string exp)
{
	string model;
	for(int i=0;i<exp.length();i++)
	{
		wchar_t ch=exp[i];
		if(ch=='('||ch==')'||ch=='+'||ch=='!')
			continue;
		if(model.find(ch)==string::npos)
			model.push_back(exp[i]);
	}

	//insertion sort
	for(int i=1;i<model.length();i++)
	{
		for(int j=i;j>0;j--)
		{
			if(model[j]<model[j-1]){
				wchar_t ch=model[j];
				model[j]=model[j-1];
				model[j-1]=ch;
			}else
				break;
		}
	}
	return model;
}

//testData is not array but a pointer, the equation of testData and model's length should be checked(added after);
bool Utility::isCNFtermFalse(string term,int testData[],string model)
{
	int result;
	
	for(int i=0;i<term.length();i++)
	{
		wchar_t ch=term[i];
		if(ch=='('||ch=='+'||ch==')')
			continue;
		if(ch=='!'){
			i++;
			ch=term[i];
			result=testData[model.find(ch)]*-1+1;
		}else
			result=testData[model.find(ch)];

		if(result==1)
			return false;
	}

	return true;
}
bool Utility::isCNFtermFalse(string term,int testData,string model)
{
	int *testArr=(int*)malloc(sizeof(int)*model.length());
	intToIntArray(testData, testArr, model.length());

	bool is=isCNFtermFalse(term,testArr,model);
	free(testArr);
	return is;
}
//testData is not array but a pointer, the equation of testData and model's length should be checked(added after);
bool Utility::isDNFtermTrue(string term,int testData[],string model)
{
	int result;
	for(int i=0;i<term.length();i++)
	{
		wchar_t ch=term[i];
		if('!'==ch){
			i++;
			ch=term[i];
			result = testData[model.find(ch)] * -1 + 1;
		}
		else{
			result = testData[model.find(ch)];
		}
		if(0==result)
			return false;
	}

	return true;
}
bool Utility::isDNFtermTrue(string term,int testData,string model)
{
	int *testArr=(int*)malloc(sizeof(int)*model.length());
	intToIntArray(testData, testArr, model.length());

	bool is=isDNFtermTrue(term,testArr,model);
	free(testArr);
	return is;
}

//exp should be DNF strictly
bool Utility::evaluateDNF(string exp,int testData[],string model)
{
	vector<string> terms;
	size_t start=0;
	size_t pos=exp.find('+');
	while(pos!=string::npos){
		terms.push_back(exp.substr(start,pos-start));
		start=pos+1;
		pos=exp.find('+',start);
	}
	terms.push_back(exp.substr(start));

	for(int i=0;i<terms.size();i++)
	{
		string term=terms.at(i);
		if(isDNFtermTrue(term,testData,model))
			return true;
	}

	return false;
}

//exp should be DNF strictly
bool Utility::evaluateDNF(string exp,int testNum, string model)
{
	int* testData=(int*)malloc(model.length()*sizeof(int));

	unsigned int base=1;
	for(int i=0;i<model.length()&&i<sizeof(int)*8;i++)
	{
		int bit=base<<i;
		if((testNum&bit)!=0)
			testData[i]=1;
		else
			testData[i]=0;
	}

	bool result=evaluateDNF(exp,testData,model);
	free(testData);
	return result;
}

//exp should be DNF strictly
bool Utility::evaluateDNF(string exp,string testString,string model)
{
	//assert(testString.length()==model.length());//de

	int* testData=(int*)malloc(model.length()*sizeof(int));
	for(int i=0;i<model.length();i++)
	{
		if(testString[i]=='0')
			testData[i]=0;
		else
			testData[i]=1;
	}

	bool result=evaluateDNF(exp,testData,model);
	free(testData);
	return result;
}

//evaluate common boolean expression
bool Utility::evaluate(string exp,int testData[],string model)
{
	if(exp.length()==0)
		return false;
	vector<string> dnfStyle=extractDNFstyle(exp);
	
	for(int i=0;i<dnfStyle.size();i++)
	{
		bool result = false;
		vector<string> cnfStyle=extractCNFstyle(dnfStyle[i]);
		for(int j=0;j<cnfStyle.size();j++)
		{
			string cnfStyleTerm=cnfStyle[j];
			if(cnfStyleTerm[0]=='!'&&cnfStyleTerm[1]=='('){
				string inBracketExp=cnfStyleTerm.substr(2,cnfStyleTerm.length()-2);
				if(inBracketExp.find('(')!=string::npos)
					result=evaluate(inBracketExp,testData,model);//recursion
				else
					result=evaluateDNF(inBracketExp,testData,model);//DNF
				result=!result;
			}else if(cnfStyleTerm[0]=='('){
				string inBracketExp=cnfStyleTerm.substr(1,cnfStyleTerm.length()-2);
				if(inBracketExp.find('(')!=string::npos)
					result=evaluate(inBracketExp,testData,model);//recursion
				else
					result=evaluateDNF(inBracketExp,testData,model);//DNF

			}else{
				result=evaluateDNF(cnfStyleTerm,testData,model);//DNF
			}
			if(!result)
				break;
		}//for
		if(result)
			return true;
	}//for

	return false;
}

string Utility::intTostring(int intNum,int length)
{
	string str;
	str.reserve(length);
	unsigned int base=1;
	for(int i=0;i<length&&i<sizeof(int)*8;i++)
	{
		int bit=base<<i;
		if((intNum&bit)!=0)
			str.push_back('1');
		else
			str.push_back('0');
	}

	return str;
}

void Utility::intToIntArray(int intNum,int intArray[],int arrayLength)
{
	unsigned int base=1;
	for(int i=0;i<arrayLength&&i<sizeof(int)*8;i++)
	{
		int bit=base<<i;
		if((intNum&bit)!=0)
			intArray[i]=1;
		else
			intArray[i]=0;
	}
}

//evaluate common boolean expression
bool Utility::evaluate(string exp,int testNum)
{
	string model=generateModel(exp);

	int* testData=(int*)malloc(model.length()*sizeof(int));

	intToIntArray(testNum,testData,model.length());

	bool result=evaluate(exp,testData,model);
	free(testData);
	return result;
}
//evaluate common boolean expression
bool Utility::evaluate(string exp, int testNum, string model)
{
	int* testData = (int*)malloc(model.length()*sizeof(int));

	intToIntArray(testNum, testData, model.length());

	bool result = evaluate(exp, testData, model);
	free(testData);
	return result;
}
//evaluate common boolean expression
bool Utility::evaluate(string exp,string testString)
{
	string model=generateModel(exp);

	//assert(testString.length()==model.length());//de

	int* testData=(int*)malloc(model.length()*sizeof(int));
	for(int i=0;i<model.length();i++)
	{
		if(testString[i]=='0')
			testData[i]=0;
		else
			testData[i]=1;
	}

	bool result=false;
	result=evaluate(exp,testData,model);
	free(testData);
	return result;
}
//evaluate mutants
bool Utility::evaluateMutant(string oriExp,string mutant,string testString)
{
	string model=generateModel(oriExp);

	//assert(testString.length()==model.length());//de

	int* testData=(int*)malloc(model.length()*sizeof(int));
	for(int i=0;i<model.length();i++)
	{
		if(testString[i]=='0')
			testData[i]=0;
		else
			testData[i]=1;
	}

	bool result=evaluate(mutant,testData,model);
	free(testData);
	return result;
}

string Utility::sortVariables(string variables)
{
	if(variables.length()==1||(variables[0]=='!'&&variables.length()==2))
		return variables;

	vector<string> var;
	for(int i=0;i<variables.length();i++)
	{
		wchar_t ch=variables[i];
		if(ch=='!'){
			i++;
			ch=variables[i];
			string s;s.push_back('!');s.push_back(ch);
			var.push_back(s);
		}else{
			string s;s.push_back(ch);
			var.push_back(s);
		}
	}//for
	
	//insertion sort
	for(int i=1;i<var.size();i++)
	{
		for(int j=i;j>0;j--)
		{
			wchar_t c1=var[j-1][0];
			if(c1=='!')
				c1=var[j-1][1];
			wchar_t c2=var[j][0];
			if(c2=='!')
				c2=var[j][1];

			if(c1>c2){
				string t=var[j];
				var[j]=var[j-1];
				var[j-1]=t;
			}else
				break;
		}//for
	}//for

	string result;
	for(int i=0;i<var.size();i++)
	{
		result.append(var[i]);
	}

	return result;
}
string Utility::sortDNF(string dnf)
{
	if(dnf.front()=='(')
		dnf=dnf.substr(1,dnf.length()-2);

	vector<string> terms=extractDNFstyle(dnf);

	for(int i=0;i<terms.size();i++)
	{
		terms[i]=sortVariables(terms[i]);
	}

	//insertion sort
	for(int i=1;i<terms.size();i++)
	{
		for(int j=i;j>0;j--)
		{
			wchar_t c1=terms[j-1][0];
			if(c1=='!')
				c1=terms[j-1][1];
			wchar_t c2=terms[j][0];
			if(c2=='!')
				c2=terms[j][1];

			if(c1>c2){
				string t=terms[j];
				terms[j]=terms[j-1];
				terms[j-1]=t;
			}else
				break;
		}//for
	}//for

	string result;
	for(int i=0;i<terms.size();i++)
	{
		result.append(terms[i]+"+");
	}
	if(result.back()=='+')
		result.pop_back();

	return result;
}
string Utility::sortCNF(string cnf)
{
	if(cnf.length()==0)
		return cnf;

	vector<string> terms=extractCNFstyle(cnf);
	string singles;
	vector<string> multiples;
	for(vector<string>::iterator itr=terms.begin();itr!=terms.end();itr++)
	{
		if(itr->length()<3)
			singles.append(*itr);
		else
			multiples.push_back(*itr);
	}

	singles=sortVariables(singles);

	//insertion sort
	for(int i=1;i<multiples.size();i++)
	{
		for(int j=i;j>0;j--)
		{
			wchar_t c1=multiples[j-1][1];
			if(c1=='!')
				c1=multiples[j-1][2];
			wchar_t c2=multiples[j][1];
			if(c2=='!')
				c2=multiples[j][2];

			if(c1>c2){
				string t=multiples[j];
				multiples[j]=multiples[j-1];
				multiples[j-1]=t;
			}else
				break;
		}//for
	}//for

	string result;
	result.append(singles);
	for(int i=0;i<multiples.size();i++)
	{
		result.append(multiples[i]);
	}

	return result;
}

bool Utility::isTwoPointsAdjacent(int p1[],int p2[],int len)
{
	int count=0;
	for(int i=0;i<len;i++)
	{
		if(p1[i]!=p2[i])
			count++;
	}

	if(count==1)
		return true;

	return false;
}

bool Utility::isTwoPointsAdjacent(int p1,int p2,int len)
{
	int* arr1=new int[len];
	int* arr2=new int[len];
	intToIntArray(p1,arr1,len);
	intToIntArray(p2,arr2,len);

	bool result=isTwoPointsAdjacent(arr1,arr2,len);

	delete[] arr1;
	delete[] arr2;

	return result;
}

vector<vector<int>> Utility::getConnectedRegion(string exp, bool trueOrFalsePoint)
{	
	Utility uti;
	string model=uti.generateModel(exp);

	vector<vector<int>> connectedRegion;
	for(int i=0;i<pow(2,model.length());i++)
	{
		if(trueOrFalsePoint&&!uti.evaluate(exp,i))
			continue;
		else if(!trueOrFalsePoint&&uti.evaluate(exp,i))
			continue;
		
		vector<int>regs;
		for(int j=0;j<connectedRegion.size();j++)
		{
			vector<int> region=connectedRegion[j];
			for(int k=0;k<region.size();k++)
			{
				if(uti.isTwoPointsAdjacent(i,region[k],model.length())){
					regs.push_back(j);
					break;
				}
			}
		}//for

		if(regs.size()==0){
			vector<int>region;
			region.push_back(i);
			connectedRegion.push_back(region);
		}else if(regs.size()==1){
			connectedRegion[regs.back()].push_back(i);
		}else{
			vector<int>&region=connectedRegion[regs.front()];
			for(int j=regs.size()-1;j>0;j--)
			{
				vector<int> second=connectedRegion[regs[j]];
				for(int k=0;k<second.size();k++)
				{
					region.push_back(second[k]);
				}
				connectedRegion.erase(connectedRegion.begin()+regs[j]);
			}//for
			region.push_back(i);
		}
	}//for
	
	vector<vector<int>> result;
	while (!connectedRegion.empty())
	{
		int p=0;
		for(int i=1;i<connectedRegion.size();i++){
			if(connectedRegion[i].size()>connectedRegion[p].size())
				p=i;
		}
		vector<int> region=connectedRegion[p];
		result.push_back(region);
		connectedRegion.erase(connectedRegion.begin()+p);
	}
	
	return result;
}


vector<int> Utility::getTermsOpCount(string exp,bool DNF1CNF0)
{
	vector<string> op;
	vector<string> up;
	vector<int> termsOpCount;

	Utility utility;
	string model=utility.generateModel(exp);
	vector<string> terms;
	if(DNF1CNF0)//
		terms=utility.extractDNFstyle(exp);//
	else
		terms=utility.extractCNFstyle(exp);//

	for (int i = 0; i < (int)pow(2, model.length()); i++) {
		int count = -1;
		for (int j = 0; j < terms.size(); j++) {
			bool is;
			if(DNF1CNF0)//
				is=utility.isDNFtermTrue(terms[j], i, model);
			else 
				is=utility.isCNFtermFalse(terms[j], i, model);
			if (is) {
				count++;
			}
		}
		if (count == 0) {
			string s=utility.intTostring(i,model.length());
			up.push_back(s);
		}else if(count>0){
			string s=utility.intTostring(i,model.length());
			op.push_back(s);
		}
	}//for

	int* a=new int[model.length()];
	for(int i=0;i<terms.size();i++)
	{
		string term=terms[i];
		for(int j=0;j<model.length();j++)
			a[j]=-1;

		for(int j=0;j<term.length();j++)
		{
			std::string::size_type pos=model.find_first_of(term[j]);
			if(pos==std::string::npos)
				continue;

			int v=1,nv=0;
			if(!DNF1CNF0){//
				v=0;nv=1;
			}
			if(j!=0&&term[j-1]=='!')
				a[pos]=nv;//
			else
				a[pos]=v;//
		}

		int opCount=0;
		for(int j=0;j<op.size();j++)
		{
			string oPoint=op[j];
			bool isCoincide=true;
			for(int k=0;k<model.length();k++)
			{
				if(a[k]==-1)
					continue;
				else if(a[k]==0&&oPoint[k]=='0')
					continue;
				else if(a[k]==1&&oPoint[k]=='1')
					continue;
				else{
					isCoincide=false;
					break;
				}
			}
			if(isCoincide)
				opCount++;
		}//for

		termsOpCount.push_back(opCount);
	
	}//for
	delete[] a;

	return termsOpCount;
}

bool Utility::isExpressionsEquivalent(string exp1, string exp2)
{
	string m1=generateModel(exp1);
	string m2=generateModel(exp2);
	if(m1.compare(m2))
		return false;

	for(int i=0;i<pow(2,m1.length());i++)
	{
		bool v1=evaluate(exp1,i);
		bool v2=evaluate(exp2,i);
		if(v1!=v2)
			return false;
	}

	return true;
}

//�Ա��ʽ����Ԥ��������ŷֵ�terms�У�vars��¼��ĸ�Ƿ����
void Utility::mutantsPreproccess(string exp,vector<string> &terms,vector<vector<string>> &terms_literals,bool (&vars)[26])
{
	Utility uti;
	terms = uti.extractDNFstyle(exp);
	for (size_t i = 0; i < terms.size(); i++)
	{
		terms_literals.push_back(uti.extractCNFstyle(terms[i]));
	}

	for (size_t i = 0; i < terms_literals.size(); i++)
	{
		for (size_t j = 0; j < terms_literals[i].size(); j++)
		{
			string temp = terms_literals[i][j];
			if (temp.size() == 1)
			{
#ifdef DEBUG
				assert(temp[0] >= 'a'&&temp[0] <= 'z');  
#endif // DEBUG

				vars[temp[0] - 'a'] = true;
			}
			else
			{
#ifdef DEBUG
				assert(temp.size() == 2 && temp[0] == '!'&& temp[1] >= 'a'&&temp[1] <= 'z');  
#endif // DEBUG

				vars[temp[1] - 'a'] = true;
			}
		}
	}
}

//�滻���ʽ�е�����,expΪԭ���ʽ��a,bΪԭ���ʽ����Ҫ���滻��term���±�
string Utility::replaceDoubleTerms(string exp, int a,int b,string term1,string term2)
{
	Utility utility;
	vector<string> terms = utility.extractDNFstyle(exp);
	int pos = exp.find(terms[a]);
	exp.replace(exp.begin()+pos,exp.begin() + pos + terms[a].size(),term1);
	pos = exp.rfind(terms[b]);
	exp.replace(exp.begin()+pos,exp.begin() + pos + terms[b].size(),term2);
	return exp;
}

//��term�ڱ��ʽ�е�λ��
int Utility::getTermPosition(vector<string> terms, string term)
{
	int pos;//term��ԭ���ʽ�����λ��
	for (int i = 0; i < terms.size(); i++){
		if (terms[i] == term){
			return i + 1;
		}
	}
	return -1;
}


//�滻���ʽ�е�����,expΪԭ���ʽ��a,bΪԭ���ʽ����Ҫ���滻��term���±�
string Utility::replaceDoubleTerms(string exp, vector<string> terms, int a, int b, string term1, string term2)
{
	Utility utility;
	int pos = exp.find(terms[a]);
	exp.replace(exp.begin() + pos, exp.begin() + pos + terms[a].size(), term1);

	pos = exp.rfind(terms[b]);
	exp.replace(exp.begin() + pos, exp.begin() + pos + terms[b].size(), term2);
	return exp;
}

bool Utility::checkTestDiffer(string test, string exprOri/*ԭ���ʽ*/, string expr/*����*/){
	Utility uti;
	bool v = (uti.evaluate(exprOri, test)) ^ (uti.evaluateMutant(exprOri, expr, test));//��򣺲�ͬΪ�棬��ͬΪ��
	return v;
}

bool Utility::initHierarchyRelation(hash_map<string, HierarchyNode>* hierarchyNodeMap,vector<HierarchyNode>* hierarchyEntry){
	
	static HierarchyNode liflif = HierarchyNode("liflif");
	static HierarchyNode liflrf = HierarchyNode("liflrf");
	static HierarchyNode liflnf = HierarchyNode("liflnf");
	//static HierarchyNode lifdorf = HierarchyNode("lifdorf");
	static HierarchyNode liftnf = HierarchyNode("liftnf");
	static HierarchyNode lifcorf = HierarchyNode("lifcorf");

	static HierarchyNode lrflrf = HierarchyNode("lrflrf");
	static HierarchyNode lrflof = HierarchyNode("lrflof");
	static HierarchyNode lrflnf = HierarchyNode("lrflnf");
	static HierarchyNode lrfcorf = HierarchyNode("lrfcorf");
	static HierarchyNode lrftnf = HierarchyNode("lrftnf");

	static HierarchyNode loflof = HierarchyNode("loflof");
	static HierarchyNode loflnf = HierarchyNode("loflnf");
	static HierarchyNode lofcorf = HierarchyNode("lofcorf");
	static HierarchyNode loftnf = HierarchyNode("loftnf");

	static HierarchyNode lnflnf = HierarchyNode("lnflnf");
	static HierarchyNode lnftnf = HierarchyNode("lnftnf");
	static HierarchyNode lnfcorf = HierarchyNode("lnfcorf");

	static HierarchyNode dorflif = HierarchyNode("dorflif");
	static HierarchyNode dorflof = HierarchyNode("dorflof");
	static HierarchyNode dorflnf = HierarchyNode("dorflnf");

	//liflif
	vector<HierarchyNode*> next;
	vector<int> nextCondition;
	next.push_back(&liflrf);
	nextCondition.push_back(1);
	//next.push_back(&lifdorf);
	//nextCondition.push_back(0);
	next.push_back(&liflnf);
	nextCondition.push_back(0);
	liflif.setNext(next);
	liflif.setNextCondition(nextCondition);
	//liflrf
	next.clear();
	nextCondition.clear();
	next.push_back(&liflnf);
	nextCondition.push_back(1);
	liflrf.setNext(next);
	liflrf.setNextCondition(nextCondition);
	//liftnf
	next.clear();
	nextCondition.clear();
	next.push_back(&liftnf);
	nextCondition.push_back(1);
	liflnf.setNext(next);
	liflnf.setNextCondition(nextCondition);


	/*

	//liflif
	vector<HierarchyNode*> next;
	vector<int> nextCondition;
	next.push_back(&liflrf);
	nextCondition.push_back(0);
	//next.push_back(&lifdorf);
	//nextCondition.push_back(0);
	next.push_back(&liflnf);
	nextCondition.push_back(0);
	liflif.setNext(next);
	liflif.setNextCondition(nextCondition);
	//liflrf
	//next.clear();
	//nextCondition.clear();
	//next.push_back(&liflnf);
	//nextCondition.push_back(2);
	//liflrf.setNext(next);
	//liflrf.setNextCondition(nextCondition);
	//liftnf
	next.clear();
	nextCondition.clear();
	next.push_back(&liftnf);
	nextCondition.push_back(0);
	liflnf.setNext(next);
	liflnf.setNextCondition(nextCondition);


	next.clear();
	nextCondition.clear();
	next.push_back(&lrflof);
	nextCondition.push_back(0);
	lrflrf.setNext(next);
	lrflrf.setNextCondition(nextCondition);

	next.clear();
	nextCondition.clear();
	next.push_back(&lrfcorf);
	nextCondition.push_back(1);
	next.push_back(&lrflnf);
	nextCondition.push_back(2);
	lrflof.setNext(next);
	lrflof.setNextCondition(nextCondition);
	




	(*hierarchyNodeMap)["liflif"] = liflif;
	(*hierarchyNodeMap)["liflrf"] = liflrf;
	(*hierarchyNodeMap)["liflnf"] = liflnf;
	(*hierarchyNodeMap)["liftnf"] = liftnf;
//	(*hierarchyNodeMap)["lifdorf"] = lifdorf;

	(*hierarchyNodeMap)["lrflrf"] = lrflrf;
	(*hierarchyNodeMap)["lrflof"] = lrflof;
	(*hierarchyNodeMap)["lrfcorf"] = lrfcorf;
	(*hierarchyNodeMap)["lrflnf"] = lrflnf;
	(*hierarchyNodeMap)["lrftnf"] = lrftnf;

	(*hierarchyEntry).push_back(liflif);
	(*hierarchyEntry).push_back(liflrf);
	(*hierarchyEntry).push_back(liflnf);
	(*hierarchyEntry).push_back(liftnf);
//	(*hierarchyEntry).push_back(lifdorf);

	(*hierarchyEntry).push_back(lrflrf);
	(*hierarchyEntry).push_back(lrflof);
	(*hierarchyEntry).push_back(lrfcorf);

	(*hierarchyEntry).push_back(lrflnf);
	*/

	(*hierarchyNodeMap)["liflrf"] = liflrf;
	(*hierarchyNodeMap)["liflnf"] = liflnf;
	(*hierarchyNodeMap)["liftnf"] = liftnf;
	return true;
}



bool Utility::getBelowNodeByCondition(int condition, string faultType, hash_map<string, HierarchyNode>& hierarchyMap, hash_map<string, HierarchyNode>& nodeSet){
	HierarchyNode thisNode = hierarchyMap[faultType];
	//��α���
	traversal(nodeSet, &thisNode, condition);
	nodeSet.erase(thisNode.getValue());
	return true;
}


bool Utility::traversal(hash_map<string, HierarchyNode>& nodeSet, HierarchyNode* node, int condition){
	if (node == nullptr){
		return true;
	}
	nodeSet[node->getValue()] == *node;
	vector<HierarchyNode *> nextList = node->getNext();
	for (int i = 0; i < nextList.size(); i++){
		if (node->getNextCondition().at(i) == condition || node->getNextCondition().at(i) == 0){
			nodeSet[nextList.at(i)->getValue()] = *nextList.at(i);
			traversal(nodeSet, nextList.at(i), condition);
		}
	}
}
