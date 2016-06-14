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
	Eigen::MatrixXf data = io->LoadData("../PSOM/src/Data/data_all_menor.csv", false);
	
	Parameter * params;
	params = new Parameter(30, 30, 300, 0.1f);
	//params = new Parameter(14, 11, 200, 0.08f);
	Som * som = new Som(data, std::make_shared<Parameter>(*params), som->RETANGULAR);

	std::cout << "The end" << std::endl;
	return 0;
}
