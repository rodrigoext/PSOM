#ifndef PARAMETER_H
#define PARAMETER_H

#include <list>
#include <string>

class Parameter
{
public:
	unsigned int map_x_, map_y_, neuron_height_, neuron_width_;
	unsigned int train_len_;
	float max_error_;
	float sigma_;
	float learning_rate_;
	float time_constant_;
	//Parameter(std::list<std::string> &params);
	Parameter(unsigned int map_size_x, unsigned int map_size_y, unsigned int train_length, float max_error);
	virtual ~Parameter();
};
#endif
