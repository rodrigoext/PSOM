#ifndef PARAMETER_H
#define PARAMETER_H
#endif

#include <list>
#include <string>

class Parameter
{
public:
	Parameter(std::list<std::string> &params);
	virtual ~Parameter();
};

