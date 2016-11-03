#include <iostream>
#include <memory>

#include "Model/Network/NeuralNetwork.h"
#include "Model/Network/Som.h"
#include "Model/Tests/Test.h"
#include "Model/Utils/IO.h"
#include "Model/Utils/Parameter.h"

#include <Eigen/Core>

int main(int argc, char *argv[])
{
	//Test * t = new Test();
	//t->TestCodebookGeneration();

	IO * io = new IO();
	Eigen::MatrixXf data = io->LoadData("../PSOM/src/Data/twod.csv", true);

	Parameter * params;
	params = new Parameter(10, 10, 15000, 0.1f, true);
	//params = new Parameter(14, 11, 200, 0.08f);
	NeuralNetwork * som = new Som(data, std::make_shared<Parameter>(*params));
	std::cout << "The end!" << std::endl;
	return 0;
}
