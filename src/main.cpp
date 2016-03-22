#include <iostream>
#include <memory>

#include "Model/Network/Som.h"
#include "Model/Tests/Test.h"
#include "Model/Utils/IO.h"
#include "Model/Utils/Parameter.h"

#include <Eigen/Dense>

int main()
{
	IO * io = new IO();
	Eigen::MatrixXf data = io->LoadData("../PSOM/src/Data/chainlink.csv");
	
	Parameter * params;
	params = new Parameter(17, 9, 100, 0.01f);
	Som * som = new Som(std::make_shared<Eigen::MatrixXf>(data), std::make_shared<Parameter>(*params));

	std::cout << "Rodou" << std::endl;
	system("PAUSE");

	return 0;
}