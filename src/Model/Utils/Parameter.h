#ifndef PARAMETER_H
#define PARAMETER_H

#include <list>
#include <string>

class Parameter
{
public:
	int map_x_, map_y_;
	int train_len_;
	float max_error_;
	//Parameter(std::list<std::string> &params);
	Parameter(int map_size_x, int map_size_y, int train_length, float max_error);
	virtual ~Parameter();
};
#endif
