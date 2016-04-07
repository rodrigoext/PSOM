#pragma once

#include "Model/Network/Som.h"

class Som;

class Train
{
public:
	enum TrainType { CLASSIC, CLASSIC_DERIVATE};
	Train(Som &som, Train::TrainType train_type = TrainType::CLASSIC);
	virtual ~Train();
private:
	bool fine_adjustment_;
	bool ClassicSomTrain(Som &som);
};

