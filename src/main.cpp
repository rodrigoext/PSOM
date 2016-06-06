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
	params = new Parameter(30, 24, 400, 0.09f);
	//params = new Parameter(14, 11, 200, 0.08f);
	Som * som = new Som(data, std::make_shared<Parameter>(*params), som->RETANGULAR);

	Eigen::MatrixXf um = som->GetUMatrix();
	io->SaveUMAT(um);

	std::cout << "The end" << std::endl;
	return 0;
}
