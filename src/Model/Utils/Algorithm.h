#pragma once
#include <vector>
#include <Eigen>

class Algorithm
{
public:
	Algorithm();
	void SetTotalEpoch(unsigned int total_epoch);
	float LearningRate(float learning_rate, unsigned int actual_epoch);
	float Radius(float sigma, unsigned int current_epoch, float time_constant);
	float CalculateNeuronDistance(Eigen::VectorXf neuron_1, Eigen::VectorXf neuron_2);
	float CalculateMedian(std::vector<float> &data);
	Eigen::MatrixXf FilterMedian(Eigen::MatrixXf &data);
	Eigen::MatrixXf Reshape(Eigen::VectorXf &data, int rows, int cols);
	float CostFunction(const Eigen::MatrixXf &data, const Eigen::MatrixXf &weigths);
	virtual ~Algorithm();
private:
	unsigned int total_epoch_;
	float get_distance(const std::vector<float> &vec1, const std::vector<float> &vec2);
};
