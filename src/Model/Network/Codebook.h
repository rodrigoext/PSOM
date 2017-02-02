#pragma once

#include <Eigen>

class Codebook
{
public:
	Codebook(unsigned int size_map_x, unsigned int size_map_y, unsigned int dimension);
	Eigen::MatrixXf GetWeights();
	void SetWeightsEndTrain(Eigen::MatrixXf new_weights);
	void Generate(Eigen::MatrixXf &data);
	virtual ~Codebook();
private:
	void Generate();

private:
	unsigned int size_;
	unsigned int dimension_;
	Eigen::MatrixXf weights_;
};
