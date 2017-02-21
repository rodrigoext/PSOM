/*
 * Watershed.h
 *
 *  Created on: 15/05/2016
 *      Author: rodrigo
 */

#pragma once

#include <vector>

#include <Eigen>

struct pos {
	int row;
	int col;
};

class Watershed {
public:
	Watershed();
	Eigen::MatrixXf transform(const Eigen::MatrixXf &input);
    Eigen::MatrixXf transform_v2(const Eigen::MatrixXf & input, int discretizarion_level = 1);
	std::vector<pos> vizinhanca4(int lin, int col, int linhas, int colunas);
	virtual ~Watershed();
};

