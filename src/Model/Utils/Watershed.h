/*
 * Watershed.h
 *
 *  Created on: 15/05/2016
 *      Author: rodrigo
 */

#pragma once

#include <vector>

#include <Eigen/Core>

struct pos {
	int l;
	int c;
};

class Watershed {
public:
	Watershed();
	Eigen::MatrixXf transform(const Eigen::MatrixXf &input);
	std::vector<pos> vizinhanca4(int lin, int col, int linhas, int colunas);
	virtual ~Watershed();
};

