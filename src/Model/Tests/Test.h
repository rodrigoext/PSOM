#ifndef TEST_H
#define TEST_H
#endif

#include "../Utils/Parameter.h"

class Test
{
public:
	Test();
	void TestLoadData();
	void TestParams(const Parameter &params);
	virtual ~Test();
};

