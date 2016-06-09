/*
 * ParetoDensity.cpp
 *
 *  Created on: 02/06/2016
 *      Author: rodrigo
 */

#include <Model/Utils/ParetoDensity.h>
#include "Algorithm.h"
#include <memory>

ParetoDensity::ParetoDensity() {
	// TODO Auto-generated constructor stub

}

ParetoDensity::~ParetoDensity() {
	// TODO Auto-generated destructor stub
}

Eigen::VectorXf ParetoDensity::CalculateDensity(Eigen::MatrixXf &data, Eigen::MatrixXf &centers, double radius)
{
	float distance;
	int dataLength = data.rows();
	Eigen::VectorXf result = Eigen::VectorXf::Zero(centers.rows());
	std::shared_ptr<Algorithm> algo;
	algo.reset(new Algorithm());
	for(int i = 0; i < centers.rows(); ++i)
	{
		auto temp = centers.row(i);
		for (int j = 0; j < dataLength; ++j)
		{
			auto tempD = data.row(j);
			float distance = (temp - tempD).squaredNorm();
			if (distance <= radius)
			{
				result(i) += 1;
			}
		}
	}
	return result;
}

