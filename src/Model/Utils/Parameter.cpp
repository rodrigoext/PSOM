#include "Parameter.h"

#include<iostream>

Parameter::Parameter(unsigned int map_size_x, unsigned int map_size_y, unsigned int train_length, float max_error)
{
	map_x_ = map_size_x;
	map_y_ = map_size_y;
	train_len_ = train_length;
	max_error_ = max_error;
}


Parameter::~Parameter()
{
}
