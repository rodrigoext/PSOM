#include "Parameter.h"

#include <iostream>
#include <algorithm>
#include <cmath>

Parameter::Parameter(unsigned int map_size_x, 
						unsigned int map_size_y, 
						unsigned int train_length, 
						float learning_rate,
						bool fine_tune)
{
	map_x_ = map_size_x;
	map_y_ = map_size_y;
	train_len_ = train_length;
	max_error_ = 0.01f;
	learning_rate_ = learning_rate;
	fine_tune_ = fine_tune;
	sigma_ = std::max(static_cast<float>(map_x_ / 2), static_cast<float>(map_y_ / 2));
	time_constant_ = train_len_ / log(sigma_);
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
