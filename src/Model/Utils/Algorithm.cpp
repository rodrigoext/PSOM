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
	return (neuron_1 - neuron_2).squaredNorm();
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

Eigen::MatrixXf Algorithm::FilterMedian(Eigen::MatrixXf &data)
{
	int rows = data.rows();
	int cols = data.cols();

	Eigen::MatrixXf result(rows, cols);
	Eigen::MatrixXf m33(3,3);
	Eigen::MatrixXf tempm33(3,3);
	Eigen::VectorXf sorted(9);

	int lv1,cv2;
	double min33;
	for (int l = 0; l < rows; l++) {
		for (int c = 0; c < cols; c++) {
			min33 = data.maxCoeff();
			for (int v1 = -1; v1 < 2; v1++)
				for (int v2 = -1; v2 < 2; v2++){
					lv1 = l+v1;
					cv2 = c+v2;
					if(lv1 >= 0 and cv2 >= 0 and lv1 < rows and cv2 < cols){
						if (data(lv1,cv2)<min33);
						min33 = data(lv1,cv2);
					}
				}

			m33.setConstant(min33);

			for (int v1 = -1; v1 < 2; v1++)
				for (int v2 = -1; v2 < 2; v2++){
					lv1 = l+v1;
					cv2 = c+v2;
					if(lv1 >= 0 and cv2 >= 0 and lv1 < rows and cv2 < cols){
						m33(v1+1,v2+1) = data(lv1,cv2);
					}
				}

			tempm33 = m33;
			double max = data.maxCoeff();
			int row_min = 0, col_min = 0;
			double temp_min = tempm33(0,0);

			for (int i = 0; i < 9; i++) {
				for (int l2 = 0; l2 < 3; l2++) {
					for (int c2 = 0; c2 < 3; c2++) {
						if (temp_min > tempm33(l2, c2)) {
							temp_min = tempm33(l2, c2);
							row_min = l2;
							col_min = c2;
						}
					}
				}
				sorted(i) = temp_min;
				tempm33(row_min, col_min) = max + 1;
				row_min = 0;
				col_min = 0;
				temp_min = tempm33(0, 0);
			}

			result(l,c) = sorted(5);
		}
	}

	return result;

};


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

float Algorithm::CostFunction(const Eigen::MatrixXf &data, const Eigen::MatrixXf &weigths)
{
	float x1, x2;
	float delta = 0;
	for (int w = 0; w < weigths.rows(); ++w)
	{
		for (int d = 0; d < data.rows(); ++d)
		{
			delta += (weigths.row(w) - data.row(d)).squaredNorm();
		}
	}
	return delta;
}

Algorithm::~Algorithm()
{
}
