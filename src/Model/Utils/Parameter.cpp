#include "Parameter.h"

#include <iostream>
#include <algorithm>

Parameter::Parameter(unsigned int map_size_x, unsigned int map_size_y, unsigned int train_length, float max_error)
{
	map_x_ = map_size_x;
	map_y_ = map_size_y;
	train_len_ = train_length;
	max_error_ = max_error;
	learning_rate_ = 0.1f;
	sigma_ = std::max(static_cast<float>(map_x_ / 2), static_cast<float>(map_y_ / 2));
	time_constant_ = train_len_ / std::log(sigma_);
}

void Parameter::SetLearningRate(float value)
{
	learning_rate_ = value;
}

void Parameter::SetSigma(float value)
{
	sigma_ = value;
}


Parameter::~Parameter()
{
}
