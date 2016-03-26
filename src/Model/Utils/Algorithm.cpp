#include "Algorithm.h"
#include <algorithm>

Algorithm::Algorithm()
{
}

void Algorithm::SetTotalEpoch(unsigned int total_epoch)
{
	total_epoch_ = total_epoch;
}

float Algorithm::LearningRate(float learning_rate, unsigned int actual_epoch)
{
	return learning_rate * std::exp(- static_cast<int>(actual_epoch) / static_cast<int>(total_epoch_));
}

float Algorithm::Radius(float sigma, unsigned int current_epoch, float time_constant)
{
	return static_cast<float>(sigma * std::exp(- static_cast<int>(current_epoch) / time_constant));
}

float Algorithm::CalculateNeuronDistance(Eigen::VectorXf neuron_1, Eigen::VectorXf neuron_2)
{
	return (neuron_1 - neuron_2).squaredNorm();
}

Algorithm::~Algorithm()
{
}
