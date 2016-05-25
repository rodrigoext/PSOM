#pragma once
#include <vector>
#include <Eigen/Core>

class Algorithm
{
public:
	Algorithm();
	void SetTotalEpoch(unsigned int total_epoch);
	float LearningRate(float learning_rate, unsigned int actual_epoch);
	float Radius(float sigma, unsigned int current_epoch, float time_constant);
	float CalculateNeuronDistance(Eigen::VectorXf neuron_1, Eigen::VectorXf neuron_2);
	float CalculateMedian(std::vector<float> &data);
	std::vector<std::vector<float>> GetUmat(int map_x, int map_y, int dimension_neurons, Eigen::MatrixXf w);
	virtual ~Algorithm();
private:
	unsigned int total_epoch_;
	float get_distance(const std::vector<float> &vec1, const std::vector<float> &vec2);
};
