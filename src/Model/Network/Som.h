#ifndef SOM_H
#define SOM_H

#include "Node.h"
#include "Model/Utils/Parameter.h"
#include "Model/Network/Codebook.h"

#include <Eigen/Core>
#include <memory>
class Som
{
private:
	std::shared_ptr<Eigen::MatrixXf> data_;
	std::shared_ptr<Parameter> params_;
	std::shared_ptr<Codebook> codebook_;
	int map_x, map_y;

public:
	Som(std::shared_ptr<Eigen::MatrixXf> data);
	Som(std::shared_ptr<Eigen::MatrixXf> data, std::shared_ptr<Parameter> params);
	virtual ~Som();
private:
	void Train();
	void CalculateMapSize();
	void DetermineRadiusInitial();
	//Node n;
};
#endif