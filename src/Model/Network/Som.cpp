#include "Som.h"
#include <iostream>
#include <limits>



Som::Som(std::shared_ptr<Eigen::MatrixXf> data)
{
	data_ = data;
	CalculateMapSize();
	//TODO: calculate map size
}

void Som::CalculateMapSize()
{
	std::cout << "Determining map size... ";
	int munits, data_len;

	if (data_->cols() <= 3)
	{
		munits = (int)(5 * pow(data_->rows(), 0.5));
		data_len = data_->rows();
		map_x = round(sqrt(munits));
		map_y = round(munits / map_x);
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


		//Eigen::EigenSolver<Eigen::MatrixXf> es(*data_.get());
		//std::cout << es.eigenvalues() << std::endl;
	}

	std::cout << "map x: " << map_x << ", map y: " << map_y << std::endl;
}

Som::~Som()
{
}
