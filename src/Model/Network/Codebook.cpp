#include "Codebook.h"

Codebook::Codebook(unsigned int size_map_x, unsigned int size_map_y, unsigned int dimension)
{
	size_ = size_map_x * size_map_y;
	dimension_ = dimension;
	Generate();
}

Eigen::MatrixXf Codebook::GetWeights()
{
	return weights_;
}

void Codebook::Generate()
{
	Eigen::MatrixXf temp(size_, dimension_);
	for (unsigned int i = 0; i < size_; ++i)
	{
		for (unsigned int j = 0; j < dimension_; ++j)
		{
			temp(i, j) = static_cast<float>(rand() / static_cast<float>(RAND_MAX));
		}
	}
	weights_ = temp;
}

Codebook::~Codebook()
{
}
