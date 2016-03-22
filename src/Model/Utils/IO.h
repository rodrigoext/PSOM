#ifndef IO_H
#define IO_H
#endif
#include <Eigen/Dense>
#include <vector>

class IO
{
public:
	IO();
	Eigen::MatrixXf LoadData(const char * file_name);
	~IO();

private:
	std::vector<std::vector<float >> ReadData(const char * file_name);
	Eigen::MatrixXf Vector2EingenMatrix(std::vector<std::vector<float>> d);
};

