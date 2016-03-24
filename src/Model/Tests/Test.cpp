#include "Test.h"

#include <iostream>

#include <Eigen/Core>
#include "Model/Network/Codebook.h"
#include "Model/Utils/IO.h"

using namespace Eigen;

Test::Test()
{
}

void Test::TestCodebookGeneration()
{
	Codebook * c = new Codebook(17, 9, 3);
	MatrixXf data = c->GetWeights();
	std::cout << data.rows() << ", " << data.cols() << std::endl;
	std::cout << data << std::endl;
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

Test::~Test()
{
}
