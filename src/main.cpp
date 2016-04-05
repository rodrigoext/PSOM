#include <iostream>
#include <memory>

#include "Model/Network/Som.h"
#include "Model/Tests/Test.h"
#include "Model/Utils/IO.h"
#include "Model/Utils/Parameter.h"

#include <Eigen/Core>

int main()
{
	//Test * t = new Test();
	//t->TestCodebookGeneration();

	IO * io = new IO();
	Eigen::MatrixXf data = io->LoadData("../PSOM/src/Data/data_seis.csv");
	
	Parameter * params;
	params = new Parameter(37, 22, 300, 0.1f);
	//params = new Parameter(14, 11, 300, 0.1f);
	Som * som = new Som(data, std::make_shared<Parameter>(*params));

	std::cout << "The end" << std::endl;
	system("PAUSE");

	return 0;
}