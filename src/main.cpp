#include <iostream>
#include <memory>

#include "Model/Network/NeuralNetwork.h"
#include "Model/Network/Som.h"
#include "Model/Tests/Test.h"
#include "Model/Utils/IO.h"
#include "Model/Utils/Parameter.h"

#include <Eigen>

int main(int argc, char *argv[])
{
	//Test * t = new Test();
	//t->TestCodebookGeneration();

	IO * io = new IO();
	Eigen::MatrixXf data = io->LoadData("src/Data/twod.csv", true);

	Parameter * params;
    params = new Parameter(15, 15, 30000, 0.1f, true);
	//params = new Parameter(14, 11, 200, 0.08f);
	Som * som = new Som(data, std::make_shared<Parameter>(*params));
	som->CalculateAllMatrix();
	std::cout << "The end!" << std::endl;
	return 0;
}
