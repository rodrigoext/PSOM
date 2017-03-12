#pragma once
#include <Eigen>
#include <vector>

class IO
{
public:
	IO();
	Eigen::MatrixXf LoadData(const char * file_name, bool normalize = true);
	void SaveUMAT(Eigen::MatrixXf &data);
	void SaveMatrix(Eigen::MatrixXf &data, std::string file_name);
	void SaveVector(Eigen::VectorXf &data, std::string file_name);
    void SaveVectorDouble(Eigen::VectorXd &data, std::string file_name);
	void MapMinMax(Eigen::MatrixXf &data);
	Eigen::MatrixXf LoadCSV(std::string file_name);
	~IO();

private:
	std::vector<std::vector<float > > ReadData(const char * file_name);
	Eigen::MatrixXf Vector2EingenMatrix(std::vector<std::vector<float> > d);
};
