/*
 * ParetoDensity.h
 *
 *  Created on: 02/06/2016
 *      Author: rodrigo
 */

#pragma once

#include <Eigen/Dense>

class ParetoDensity {
public:
	ParetoDensity();
	Eigen::VectorXf CalculateDensity(Eigen::MatrixXf &data, Eigen::MatrixXf &centers, double radius);
	virtual ~ParetoDensity();
};
