#pragma once
#include <Eigen/Dense>
#include <vector>

class IO
{
public:
	IO();
	Eigen::MatrixXf LoadData(const char * file_name);
	~IO();

private:
	std::vector<std::vector<float > > ReadData(const char * file_name);
	void MapMinMax(Eigen::MatrixXf &data);
	Eigen::MatrixXf Vector2EingenMatrix(std::vector<std::vector<float> > d);
};

