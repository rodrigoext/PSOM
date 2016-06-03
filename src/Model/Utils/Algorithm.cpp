#include "Algorithm.h"
#include <algorithm>
#include <iostream>

using namespace std;

Algorithm::Algorithm()
{
}

void Algorithm::SetTotalEpoch(unsigned int total_epoch)
{
	total_epoch_ = total_epoch;
}

float Algorithm::LearningRate(float learning_rate, unsigned int actual_epoch)
{
	return static_cast<float>(learning_rate * std::exp(- static_cast<int>(actual_epoch) / static_cast<int>(total_epoch_)));
}

float Algorithm::Radius(float sigma, unsigned int current_epoch, float time_constant)
{
	return static_cast<float>(sigma * std::exp(- static_cast<int>(current_epoch) / time_constant));
}

float Algorithm::CalculateNeuronDistance(Eigen::VectorXf neuron_1, Eigen::VectorXf neuron_2)
{
	neuron_1.transpose();
	neuron_2.transpose();
	float result = (neuron_1 - neuron_2).dot(neuron_1 - neuron_2);
	return sqrt(result);
}

float Algorithm::CalculateMedian(std::vector<float> &data)
{
	float median = 0.0f;
	size_t size = data.size();
	sort(data.begin(), data.end());
	if (size % 2 == 0)
	{
		median = (data[size / 2 - 1] + data[size / 2]) / 2;
	}
	else
	{
		median = data[size / 2];
	}
	return median;
}

float Algorithm::get_distance(const std::vector<float> &vec1, const std::vector<float> &vec2)
{
	float distance = 0;

	for (int i = 0; i < vec2.size(); i++)
	{
		distance += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
	}

	return sqrt(distance);
	//return distance;
}

Eigen::MatrixXf Algorithm::Reshape(Eigen::VectorXf &data, int rows, int cols)
{
	return Eigen::Map<Eigen::MatrixXf>(data.data(), rows, cols);
}

Algorithm::~Algorithm()
{
}
