#pragma once
#include <Eigen/Core>

class Algorithm
{
public:
	Algorithm();
	void SetTotalEpoch(unsigned int total_epoch);
	float LearningRate(float learning_rate, unsigned int actual_epoch);
	float Radius(float sigma, unsigned int current_epoch, float time_constant);
	float CalculateNeuronDistance(Eigen::VectorXf neuron_1, Eigen::VectorXf neuron_2);
	virtual ~Algorithm();
private:
	unsigned int total_epoch_;
};
