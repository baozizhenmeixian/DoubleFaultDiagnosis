#pragma once

#include <string>
#include <vector>
#include <hash_set>
#include <hash_map>
#include "Mutants.h"
#include "HierarchyNode.h"

class DnfMutant
{
public:
	DnfMutant(void);
	~DnfMutant(void);

	void addMutant(string exp, string type);

	Mutants getMutants(string exp, string faultType);
	Mutants getMutants(string exp, string faultType, vector<string> terms);

	Mutants getMutants(string exp, string faultType, vector<string> oriTrems, vector<string> terms, string model);

	vector<string> subTermOver2Exp(string exp, vector<string>terms, int start, int end);


	//============================================== 双缺陷——	按项生成	==============================================
	//【一】词 —— 19个，全部无序
	Mutants generateLNFxLNFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLNF---(1)
	Mutants generateDLNNFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLNF---(2)

	Mutants generateLNFxLOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLOF---(3)
	Mutants generateDLNOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLOF---(4)

	Mutants generateLNFxLIFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLIF---(5)
	Mutants generateDLNIFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLIF---(6)

	Mutants generateLNFxLRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLRF---(7)
	Mutants generateDLNRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LNFxLRF---(8)

	Mutants generateLOFxLOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LOFxLOF---(9)
	Mutants generateDLOOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LOFxLOF---(10)

	Mutants generateLOFxLIFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LOFxLIF---(11)

	Mutants generateLOFxLRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LOFxLRF---(12)
	Mutants generateDLORFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LOFxLRF ---(13)

	Mutants generateLIFxLIFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LIFxLIF---(14)
	Mutants generateDLIIFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LIFxLIF---(15)

	Mutants generateLIFxLRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LIFxLRF---(16)
	Mutants generateDLIRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LIFxLRF---(17)

	Mutants generateLRFxLRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //LRFxLRF---(18)
	Mutants generateDLRRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//LRFxLRF---(19)


	//【二】项 —— 27个无序、4个有序
	// 无序
	//ENF相关----（1）-（12）
	Mutants generateENFxENFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxENF---(1)

	Mutants generateENFxTNFCase1doublemutants(string subexp1, string subexp2, string faultType, string model); //ENFxTNF---(3)
	Mutants generateENFxTNFCase2doublemutants(string exp, int i1, int h1, int i2, int h2, string faultType, string model); //ENFxTNF---(4)

	Mutants generateENFxTNFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxTNF---(3)

	Mutants generateENFxTOFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxTOF---(5)

	Mutants generateENFxDORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxDORF---(7)

	Mutants generateENFxCORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxCORF---(11)

	//-------------------------------------------------------------------------------------------------------------------------------
	Mutants generateTNFxTNFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxTNF---(13)
	Mutants generateTNFxTOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxTOF---(14)

	Mutants generateTNFxDORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//TNFxDORF---(15)—两个相邻的项发生DORF，剩下的项发生TNF
	Mutants generateTNFxDORFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//TNFxDORF---(16)—term1，term2中有一个取反，然后合并（term1必须与term2相邻）

	Mutants generateTNFxCORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//TNFxCORF---(17)
	Mutants generateTNFxCORFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model);//TNFxCORF---(18)

	Mutants generateTOFxTOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxTOF---(19)

	Mutants generateTOFxDORFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxDORF---(20)
	Mutants generateTOFxCORFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxCORF---(21)

	Mutants generateDORFxDORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxDORF---(22)
	Mutants generateDORFxDORFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxDORF---(23)

	Mutants generateDORFxCORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxCORF---(24)——前两个项DORF,另一个项CORF
	Mutants generateDORFxCORFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxCORF---(25)——前一个项CORF,分开的后一半与另一个项DORF

	Mutants generateCORFxCORFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxCORF---(26)
	Mutants generateCORFxCORFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxCORF---(27)

	// 有序
	Mutants generateTOFxDORFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxDORF---(53) term2遗漏，term1和term3合并
	Mutants generateCORFxTNFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxTNF---(73)
	Mutants generateCORFxTOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxTOF---(76)

	//【三】 词 + 项 ~~~~~~~~~~~ 
	Mutants generateENFxLNFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxLNF---(1)

	Mutants generateENFxLOFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxLOF---(3)

	Mutants generateENFxLIFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxLIF---(5)

	Mutants generateENFxLRFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //ENFxLIF---(7)

	Mutants generateTNFxLNFTwoTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLNF---(9)
	Mutants generateTNFxLNFOneTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLNF---(10)

	Mutants generateTNFxLOFTwoTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLOF---(11)
	Mutants generateTNFxLOFOneTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLOF---(12)

	Mutants generateTNFxLIFTwoTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLIF---(13)
	Mutants generateTNFxLIFOneTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLIF---(14)

	Mutants generateTNFxLRFTwoTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLRF---(15)
	Mutants generateTNFxLRFOneTermdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLRF---(16)

	Mutants generateTOFxLNFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxLNF---(17)
	Mutants generateTOFxLOFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxLOF---(18)
	Mutants generateTOFxLIFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxLIF---(19)
	Mutants generateTOFxLRFdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TOFxLRF---(20)

	Mutants generateDORFxLNFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxLNF---(21)
	Mutants generateDORFxLNFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxLNF---(22)

	Mutants generateDORFxLOFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxLOF---(23)
	Mutants generateDORFxLOFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> faultterms, string model); //DORFxLOF---(24)

	Mutants generateDORFxLIFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxLIF---(25)
	Mutants generateDORFxLIFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> faultterms, string model); //DORFxLIF---(26)

	Mutants generateDORFxLRFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //DORFxLRF---(27)
	Mutants generateDORFxLRFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> faultterms, string model); //DORFxLRF---(28)

	Mutants generateCORFxLNFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLNF---(29)
	Mutants generateCORFxLNFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLNF---(30)

	Mutants generateCORFxLOFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLOF---(31)
	Mutants generateCORFxLOFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLOF---(32)

	Mutants generateCORFxLIFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLIF---(33)
	Mutants generateCORFxLIFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLIF---(34)

	Mutants generateCORFxLRFCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLRF---(35)
	Mutants generateCORFxLRFCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLRF---(36)

	//有序
	Mutants generateTNFxLIFWithOrderdoublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //TNFxLIF---(52)

	Mutants generateCORFxLIFWithOrderCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLIFWithOrderCase1---(73)
	Mutants generateCORFxLIFWithOrderCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLIFWithOrderCase2---(74)

	Mutants generateCORFxLRFWithOrderCase1doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLRFWithOrderCase1---(77)
	Mutants generateCORFxLRFWithOrderCase2doublemutants(string exp, vector<string> oriTerms, vector<string> terms, string model); //CORFxLRFWithOrderCase2---(78)


	//==============================================================================================


	//单缺陷
	Mutants generateLIFmutants(string exp);
	Mutants generateLOFmutants(string exp);
	Mutants generateLRFmutants(string exp);
	Mutants generateLNFmutants(string exp);
	Mutants generateTOFmutants(string exp);
	Mutants generateTNFmutants(string exp);
	Mutants generateCORFmutants(string exp);
	Mutants generateDORFmutants(string exp);
	Mutants generateENFmutants(string exp);
	Mutants generateLIORFmutants(string exp);

	Mutants generateLIFmutants(string exp, vector<string> terms);
	Mutants generateLOFmutants(string exp, vector<string> terms);
	Mutants generateLNFmutants(string exp, vector<string> terms);
	Mutants generateLRFmutants(string exp, vector<string> terms);
	Mutants generateTOFmutants(string exp, vector<string> terms);
	Mutants generateTNFmutants(string exp, vector<string> terms);
	Mutants generateCORFmutants(string exp, vector<string> terms);
	Mutants generateDORFmutants(string exp, vector<string> terms);
	Mutants generateLIORFmutants(string exp, vector<string> terms);


	//=== 单项双缺陷 ===
	//词相关：
	Mutants generateDLIRFdoublemutants(string exp);//LIF x LRF --
	Mutants generateDLRRFdoublemutants(string exp);//LRF x LRF --
	Mutants generateDLORFdoublemutants(string exp);//LOF x LRF --

	Mutants generateSingleTermLIFxLIFdoublemutants(string exp);//LIF x LIF
	Mutants generateSingleTermLOFxLOFdoublemutants(string exp);//LOF x LOF
	Mutants generateSingleTermLNFxLNFdoublemutants(string exp);//LNF x LNF
	Mutants generateSingleTermLIFxLNFdoublemutants(string exp);//LIF x LNF
	//LIFxLOF等价于LRF
	Mutants generateSingleTermLOFxLNFdoublemutants(string exp);//LOF x LNF
	Mutants generateSingleTermLNFxLRFdoublemutants(string exp);//LNF x LRF

	//词+项(TNF)相关：
	Mutants generateSingleTermLIFxTNFdoublemutants(string exp);//LIF x TNF
	Mutants generateSingleTermLRFxTNFdoublemutants(string exp);//LRF x TNF 
	Mutants generateSingleTermLOFxTNFdoublemutants(string exp);//LOF x TNF
	Mutants generateSingleTermLNFxTNFdoublemutants(string exp);//LNF x TNF

	//CORF + 词 相关：
	Mutants generateSingleTermCORFxLIFmutants(string exp);
	Mutants generateSingleTermCORFxLOFmutants(string exp);
	Mutants generateSingleTermCORFxLRFmutants(string exp);
	Mutants generateSingleTermCORFxLNFmutants(string exp);

	//CORF+TNF

	//DORF + 词 相关：
	Mutants generateSingleTermDORFxLIFmutants(string exp);
	Mutants generateSingleTermDORFxLOFmutants(string exp);
	Mutants generateSingleTermDORFxLNFmutants(string exp);



	//=== 双项双缺陷 ===
	Mutants generateLIFxLIFdoublemutants(string exp);//LIF x LIF
	Mutants generateLIFxLRFdoublemutants(string exp);//LIF x LRF
	Mutants generateLRFxLRFdoublemutants(string exp);//LRF x LRF
	Mutants generateDORFxLRFdoublemutants(string exp);//DORF x LRF
	Mutants generateTOFxLRFdoublemutants(string exp);//TOFxLRF

	Mutants DnfMutant::diagnosisSingleTermLIFxLIFdoublemutants(
		string exp,
		string faultExp,//待测表达式
		hash_set<string>& outMutant,//排除的变体格式
		vector<vector<string>> optiUniformitySet,//最优无效点
		vector<vector<string>> optiDifferSet,//最优有效点
		hash_map<string, HierarchyNode>& hierarchyMap
		);

	Mutants DnfMutant::diagnosisSingleTermLIFxLRFdoublemutants(
		string exp,
		string faultExp,//待测表达式
		hash_set<string>& outMutant,//排除的变体格式
		vector<vector<string>> optiUniformitySet,//最优无效点
		vector<vector<string>> optiDifferSet,//最优有效点
		hash_map<string, HierarchyNode>& hierarchyMap
		);


	//对表达式进行预处理，将项才分到terms中，vars记录字母是否出现
	//void mutantsPreproccess(string exp, vector<string> &terms, bool(&vars)[26]);
};
