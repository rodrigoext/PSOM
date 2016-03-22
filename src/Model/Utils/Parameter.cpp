#include "Parameter.h"

#include<iostream>

Parameter::Parameter(std::list<std::string> &params)
{
	for (std::list<std::string>::iterator it = params.begin(); it != params.end(); ++it)
	{
		std::cout << *it << std::endl;
	}
}


Parameter::~Parameter()
{
}
