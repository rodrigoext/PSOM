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

	std::cout << neuron_1.transpose() << std::endl;
	std::cout << neuron_2.transpose() << std::endl;
	neuron_1.transpose();
	neuron_2.transpose();
	float result = (neuron_1 - neuron_2).dot(neuron_1 - neuron_2);
	std::cout << "--" << std::endl;
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

std::vector<std::vector<float>> Algorithm::GetUmat(int map_x, int map_y, int dimension_neurons, Eigen::MatrixXf w) {
	vector<float> um;
	vector<vector<vector<float>>> neurons(map_x, vector<vector<float>>(map_y, vector<float>(dimension_neurons)));

	int count = 0;

	/*for (int i = 0; i < map_x; ++i)
	{
		for (int j = 0; j < map_y; j++)
		{
			vector<float> temp;
			for (int k = 0; k < dimension_neurons; ++k)

			neurons[i][j] = w.row(count);
			count++;
		}
	}*/

	int X = 2 * map_x - 1;
	int Y = 2 * map_y - 1;

	vector<vector<float>> U(X, vector<float>(Y));

	for (int i = 0; i < map_x; ++i)
	{
		for (int j = 0; j < map_y; ++j)
		{
			if (j < (map_y - 1))
			{
				float dx = get_distance(neurons[i][j], neurons[i][j + 1]);
				U[2 * i][2 * j + 1] = dx;
				um.push_back(dx);
			}

			if (i < (map_x - 1))
			{
				float dy = get_distance(neurons[i][j], neurons[i + 1][j]);
				U[2 * i + 1][2 * j] = dy;
				um.push_back(dy);
			}

			if (i < (map_x - 1) && j < (map_y - 1))
			{
				float dz1 = get_distance(neurons[i][j], neurons[i+1][j+1]);
				float dz2 = get_distance(neurons[i+1][j], neurons[i][j+1]);
				U[2 * i + 1][2 * j + 1] = (dz1 + dz2) / (2 * sqrt(2));
			}
		}

		vector<float> a;

		for (int i = 0; i < X; i += 2)
		{
			for (int j = 0; j < Y; j += 2)
			{
				if (j > 0 && i > 0 && i < (X - 1) && j < (Y - 1))
				{
					a.push_back(U[i][j-1]);
					a.push_back(U[i][j+1]);
					a.push_back(U[i-1][j]);
					a.push_back(U[i+1][j]);
				}
				else if (i == 0 && j > 0 && j < (Y - 1))
				{
					a.push_back(U[i][j-1]);
					a.push_back(U[i][j+1]);
					a.push_back(U[i+1][j]);
				}
				else if (i == (X - 1) && j > 0 && j < (Y - 1))
				{
					a.push_back(U[i][j - 1]);
					a.push_back(U[i][j + 1]);
					a.push_back(U[i - 1][j]);
				}
				else if (j == 0 && i > 0 && i < (X - 1))
				{
					a.push_back(U[i][j + 1]);
					a.push_back(U[i - 1][j]);
					a.push_back(U[i + 1][j]);
				}
				else if (j == (Y - 1) && i > 0 && i < (X - 1))
				{
					a.push_back(U[i][j - 1]);
					a.push_back(U[i - 1][j]);
					a.push_back(U[i + 1][j]);
				}
				else if (j == 0 && i == 0)
				{
					a.push_back(U[i][j + 1]);
					a.push_back(U[i + 1][j]);
				}
				else if (i == 0 && j > 0 && j < (Y - 1))
				{
					a.push_back(U[i][j]);
					a.push_back(U[i][j]);
					a.push_back(U[i][j]);
				}
				else if (j == (X - 1) && i == 0)
				{
					a.push_back(U[i][j - 1]);
					a.push_back(U[i + 1][j]);
				}
				else if (j == 0 && i < (X - 1))
				{
					a.push_back(U[i][j + 1]);
					a.push_back(U[i - 1][j]);
				}
				else if (j == (Y - 1) && i == (X - 1))
				{
					a.push_back(U[i][j - 1]);
					a.push_back(U[i - 1][j]);
				}
				else
				{
					a.push_back(0.0f);
				}

				U[i][j] = CalculateMedian(a);
			}
		}
	}

	return U;
}

Algorithm::~Algorithm()
{
}
