#pragma once
#include"math.h"
class CCell
{
public:
	CCell(void);
	~CCell(void);
	CCell(int dimension);
	void Copy(CCell* pCell);
	int dimension;
	int* value;//dimension�����������������
	//���cell��tModel�е�������311
	int getTModelIndex();
	bool compare(CCell* pCell);
	
};

