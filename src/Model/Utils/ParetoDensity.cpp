/*
 * ParetoDensity.cpp
 *
 *  Created on: 02/06/2016
 *      Author: rodrigo
 */

#include <Model/Utils/ParetoDensity.h>
#include "Algorithm.h"
#include <memory>
#include <iostream>

ParetoDensity::ParetoDensity() {
	// TODO Auto-generated constructor stub

}

ParetoDensity::~ParetoDensity() {
	// TODO Auto-generated destructor stub
}

Eigen::VectorXd ParetoDensity::CalculateDensity(Eigen::MatrixXf &data, Eigen::MatrixXf &centers, double radius)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(centers.rows());
	std::shared_ptr<Algorithm> algo;
	algo.reset(new Algorithm());
    #pragma omp parallel for
    for(int i = 0; i < centers.rows(); ++i) {
        for (int j = 0; j < data.rows(); ++j) {
            if ((centers.row(i) - data.row(j)).squaredNorm() <= radius) {
				result(i) += 1;
			}
		}
	}
	return result;
}

