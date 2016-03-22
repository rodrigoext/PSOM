#include "Test.h"

#include <iostream>

#include "../Utils/IO.h"
#include "Eigen/Dense"


using namespace Eigen;

Test::Test()
{
}

void Test::TestLoadData()
{
	IO * io;
	io = new IO();
	MatrixXf data = io->LoadData("../PSOM/src/Data/chainlink.csv");

	data.normalize();

	std::cout << data << std::endl;

	delete io;
}

void Test::TestParams(const Parameter &params)
{

}

Test::~Test()
{
}
