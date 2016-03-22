#include <iostream>
#include <memory>

#include "Model/Network/Som.h"
#include "Model/Tests/Test.h"
#include "Model/Utils/IO.h"

#include <Eigen/Dense>

int main()
{
	Test * t;
	t = new Test();

	//t->TestLoadData();
	IO * io = new IO();
	Eigen::MatrixXf data = io->LoadData("../PSOM/src/Data/chain.csv");
	
	Som * som = new Som(std::make_shared<Eigen::MatrixXf>(data));

	std::cout << "Rodou" << std::endl;
	system("PAUSE");

	return 0;
}