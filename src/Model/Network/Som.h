#pragma once

#include "Node.h"
#include "Model/Utils/Parameter.h"
#include "Model/Utils/Algorithm.h"
#include "Model/Network/Codebook.h"
#include "Model/Network/Train.h"

#include <Eigen/Core>
#include <memory>
class Som
{
	friend class Train;
protected:
	Eigen::MatrixXf data_;
	Eigen::MatrixXf grid_;
	Eigen::MatrixXf umat_;
	std::shared_ptr<Parameter> params_;
	std::shared_ptr<Codebook> codebook_;
	std::shared_ptr<Algorithm> algorithm_;
	int map_x, map_y;

public:
	enum Topology
	{
		RETANGULAR,
		HEXAGONAL
	};
	Som(Eigen::MatrixXf data, Som::Topology topology = Topology::HEXAGONAL);
	Som(Eigen::MatrixXf data, std::shared_ptr<Parameter> params, Som::Topology topology = Topology::HEXAGONAL);
	virtual ~Som();
private:
	void InitGrid(Som::Topology topology);
	void TrainSom();
	void CalculateMapSize();
	void DetermineRadiusInitial();
	void NInv(int n, int &width, int &height);
	void CalculateUMatrix();
	//Node n;
};