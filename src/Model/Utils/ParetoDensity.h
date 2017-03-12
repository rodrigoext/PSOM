/*
 * ParetoDensity.h
 *
 *  Created on: 02/06/2016
 *      Author: rodrigo
 */

#pragma once

#include <Eigen>

class ParetoDensity {
public:
	ParetoDensity();

    Eigen::VectorXd CalculateDensity(Eigen::MatrixXf &data, Eigen::MatrixXf &centers, double radius);
	virtual ~ParetoDensity();
};
