#include "Codebook.h"
#include <random>
#include <iostream>
#include <Model/Utils/IO.h>
#include <Model/Utils/Algorithm.h>

Codebook::Codebook(unsigned int size_map_x, unsigned int size_map_y, unsigned int dimension)
{
	size_ = size_map_x * size_map_y;
	dimension_ = dimension;
	//Generate();
}

Eigen::MatrixXf Codebook::GetWeights()
{
	return weights_;
}

void Codebook::SetWeightsEndTrain(Eigen::MatrixXf new_weights)
{
	weights_ = new_weights;
}

void Codebook::Generate()
{
	Eigen::MatrixXf temp(size_, dimension_);
	std::random_device rd;
	std::mt19937 eng(rd());
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	for (unsigned int i = 0; i < size_; ++i)
	{
		for (unsigned int j = 0; j < dimension_; ++j)
		{
			//temp(i, j) = static_cast<float>(rand() / static_cast<float>(RAND_MAX));
			temp(i, j) = static_cast<float>(dis(eng));
			//temp(i, j) = 0.0f;
		}
	}
	IO * io = new IO();
	io->SaveMatrix(temp, "codebook_init");
	delete io;
	weights_ = temp;
}

void Codebook::Generate(Eigen::MatrixXf &data)
{
	Eigen::MatrixXf temp(size_, dimension_);
	Eigen::MatrixXf temp_best(size_, dimension_);
	std::random_device rd;
	std::mt19937 eng(rd());
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	Algorithm *algo = new Algorithm();
	float best_delta = 9999999.0f;
	float temp_delta;
	for (int k = 0; k < 10; ++k)
	{
		for (unsigned int i = 0; i < size_; ++i)
			for (unsigned int j = 0; j < dimension_; ++j)
				temp(i, j) = static_cast<float>(dis(eng));

		temp_delta = algo->CostFunction(data, temp);
		if (temp_delta < best_delta)
		{
			best_delta = temp_delta;
			temp_best = temp;
		}
	}
	weights_ = temp_best;
	IO * io = new IO();
	io->SaveMatrix(temp_best, "codebook_init");
	delete io;
}

Codebook::~Codebook()
{
}
