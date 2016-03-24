#include "Som.h"
#include <iostream>
#include <limits>

Som::Som(std::shared_ptr<Eigen::MatrixXf> data)
{
	data_ = data;
	CalculateMapSize();
	DetermineRadiusInitial();
	//TODO: calculate map size
}

Som::Som(std::shared_ptr<Eigen::MatrixXf> data, std::shared_ptr<Parameter> params)
{
	data_ = data;
	params_ = params;
	codebook_.reset(new Codebook(params_->map_x_, params_->map_y_, data_->cols()));
	std::cout << codebook_->GetWeights() << std::endl;
}

void Som::CalculateMapSize()
{
	std::cout << "Determining map size... ";
	int munits, data_len;

	if (data_->cols() <= 3)
	{
		munits = (int)(5 * pow(data_->rows(), 0.5));
		data_len = data_->rows();
		map_x = static_cast<int>(round(sqrt(munits)));
		map_y = static_cast<int>(round(munits / map_x));
	}
	else
	{
		Eigen::MatrixXf A(data_->cols(), data_->cols());
		A.setZero();
		for (int i = 0; i < A.rows(); ++i)
			for (int j = 0; j < A.cols(); ++j)
				A(i,j) = std::numeric_limits<float>::infinity();

		for (int i = 0; i < A.rows(); ++i)
		{
			for (int j = i; j < A.cols(); ++j)
			{
				Eigen::VectorXf op1 = data_->col(i);
				std::cout << op1 << std::endl;
				Eigen::VectorXf op2 = data_->col(j);
				Eigen::MatrixXf op3 = op1*op2.transpose();
				std::cout << op3<< std::endl;
				for (int k = 0; k < data_->rows(); ++k)
				{

				}
			}
		}
	}
	std::cout << "map x: " << map_x << ", map y: " << map_y << std::endl;
}

Som::~Som()
{
}

void Som::DetermineRadiusInitial()
{
	int ms = std::max(map_x, map_y);
	float map_radius_ini = fmaxf(1.0f, ((float)ms / 4));
	float map_radius_fin = fmaxf(1, (map_radius_ini / 4));
	int train_len = std::max(1, (int)(10 * ((map_x * map_y) / data_->rows())));
}
