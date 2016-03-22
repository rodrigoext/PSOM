#ifndef SOM_H
#define SOM_H
#endif

#include "Node.h"
//#include "../Model/Utils/Parameter.h"

#include <Eigen/Dense>
#include <memory>
class Som
{
private:
	std::shared_ptr<Eigen::MatrixXf> data_;
	int map_x, map_y;

public:
	Som(std::shared_ptr<Eigen::MatrixXf> data);
	//Som(Eigen::MatrixXf *data, Parameter *params);
	virtual ~Som();
private:
	void CalculateMapSize();
	//Node n;
};
