#ifndef SOM_H
#define SOM_H

#include "Node.h"
#include "Model/Utils/Parameter.h"
#include "Model/Utils/Algorithm.h"
#include "Model/Network/Codebook.h"

#include <Eigen/Core>
#include <memory>
class Som
{
private:
	Eigen::MatrixXf data_;
	Eigen::MatrixXf grid_;
	std::shared_ptr<Parameter> params_;
	std::shared_ptr<Codebook> codebook_;
	std::shared_ptr<Algorithm> algorithm_;
	int map_x, map_y;

public:
	Som(Eigen::MatrixXf data);
	Som(Eigen::MatrixXf data, std::shared_ptr<Parameter> params);
	virtual ~Som();
private:
	void InitGrid();
	void Train();
	void CalculateMapSize();
	void DetermineRadiusInitial();
	void NInv(int n, int &width, int &height);
	//Node n;
};
#endif