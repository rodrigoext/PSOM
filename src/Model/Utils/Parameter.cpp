#include "Parameter.h"

#include<iostream>

Parameter::Parameter(int map_size_x, int map_size_y, int train_length, float max_error)
{
	map_x_ = map_size_x;
	map_y_ = map_size_y;
	train_len_ = train_length;
	max_error_ = max_error_;
}


Parameter::~Parameter()
{
}
