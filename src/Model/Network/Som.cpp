#include "Som.h"
#include <iostream>
#include <limits>
#include <vector>

#include <Model/Utils/Watershed.h>
#include <Model/Utils/ParetoDensity.h>
#include <Model/Utils/IO.h>

Som::Som(Eigen::MatrixXf data, Som::Topology topology)
{
	data_ = data;
	CalculateMapSize();
	DetermineRadiusInitial();
	//TODO: calculate map size
}

Som::Som(Eigen::MatrixXf data, std::shared_ptr<Parameter> params, Som::Topology topology)
{
	data_ = data;
	params_ = params;
	map_x = params_->map_x_;
	map_y = params_->map_y_;
	algorithm_.reset(new Algorithm());
	algorithm_->SetTotalEpoch(params_->train_len_);
	codebook_.reset(new Codebook(params_->map_x_, params_->map_y_, data_.cols()));
	//std::cout << codebook_->GetWeights() << std::endl;
	//std::cout << "---------------------------------------" << std::endl;
	InitGrid(topology);
	std::shared_ptr<Train> t(new Train(*this));
	CalculateUMatrix();
	CalculatePMatrix();
	//TrainSom();
}

void Som::InitGrid(Som::Topology topology)
{
	grid_ = Eigen::MatrixXf(codebook_->GetWeights().rows(), 2);
	int width, height;
	switch (topology)
	{
	case HEXAGONAL:
		for (int i = 1; i <= map_x*map_y; i++) 
		{
			NInv(i, height, width);
			grid_(i - 1, 0) = static_cast<float>(width);
			grid_(i - 1, 1) = static_cast<float>(height) * 0.8660;
			if (height % 2 == 0)
				grid_(i - 1, 0) += 0.5;
		}
		break;
	case RETANGULAR:
		for (int i = 1; i <= map_x*map_y; ++i)
		{
			NInv(i, height, width);
			grid_(i - 1, 0) = static_cast<float>(width);
			grid_(i - 1, 1) = static_cast<float>(height);
		}
		break;
	default:
		break;
	}
	std::cout << "grid -------" << std::endl;
	std::cout << grid_ << std::endl;
	std::cout << "fim grid -------" << std::endl;
}

void Som::TrainSom()
{
	std::cout << "Training... " << std::endl;
	float sigma = params_->sigma_;
	//Eigen::VectorXf numerator(map_x * map_y * data_.cols());
	//Eigen::VectorXf denominator(map_x * map_y);
	int dim = data_.cols();
	float dist, learning_rate = 0.0f;
	float hci_exp = 0.0f;
	Eigen::MatrixXf::Index index;
	Eigen::MatrixXf weights = codebook_->GetWeights();
	for (int current_epoch = 0; current_epoch < params_->train_len_; ++current_epoch)
	{
		sigma = algorithm_->Radius(params_->sigma_, current_epoch, params_->time_constant_);
		learning_rate = algorithm_->LearningRate(params_->learning_rate_, current_epoch);
		int rand_sample = rand() % (data_.rows() + 1);
		auto sample = data_.row(rand_sample);
		//get BMU
		(weights.rowwise() - sample).rowwise().squaredNorm().minCoeff(&index);
		auto bmu = grid_.row(index);
		for (int n = 0; n < weights.rows(); ++n)
		{
			dist = (bmu - grid_.row(n)).squaredNorm();
			float sigma2 = sigma*sigma;
			if (dist < sigma2)
			{
				hci_exp = std::exp((- dist) / (2 * sigma2));
				weights.row(n) += learning_rate * hci_exp * (sample - weights.row(n));
			}
			
		}
	}
	//Ajuste fino
	params_->SetSigma(0.1f);
	params_->SetLearningRate(0.01f);

	for (int current_epoch = 0; current_epoch < params_->train_len_; ++current_epoch)
	{
		sigma = algorithm_->Radius(params_->sigma_, current_epoch, params_->time_constant_);
		learning_rate = algorithm_->LearningRate(params_->learning_rate_, current_epoch);
		int rand_sample = rand() % (data_.rows() + 1);
		auto sample = data_.row(rand_sample);
		//get BMU
		(weights.rowwise() - sample).rowwise().squaredNorm().minCoeff(&index);
		auto bmu = grid_.row(index);
		for (int n = 0; n < weights.rows(); ++n)
		{
			dist = (bmu - grid_.row(n)).squaredNorm();
			float sigma2 = sigma*sigma;
			if (dist < sigma2)
			{
				hci_exp = std::exp((-dist) / (2 * sigma2));
				weights.row(n) += learning_rate * hci_exp * (sample - weights.row(n));
			}

		}
	}

	codebook_->SetWeightsEndTrain(weights);
	std::cout << "Final weights" << std::endl;
	std::cout << codebook_->GetWeights() << std::endl;
}

void Som::CalculateMapSize()
{
	std::cout << "Determining map size... ";
	int munits, data_len;

	if (data_.cols() <= 3)
	{
		munits = (int)(5 * pow(data_.rows(), 0.5));
		data_len = data_.rows();
		map_x = static_cast<int>(round(sqrt(munits)));
		map_y = static_cast<int>(round(munits / map_x));
	}
	else
	{
		Eigen::MatrixXf A(data_.cols(), data_.cols());
		A.setZero();
		for (int i = 0; i < A.rows(); ++i)
			for (int j = 0; j < A.cols(); ++j)
				A(i,j) = std::numeric_limits<float>::infinity();

		for (int i = 0; i < A.rows(); ++i)
		{
			for (int j = i; j < A.cols(); ++j)
			{
				Eigen::VectorXf op1 = data_.col(i);
				std::cout << op1 << std::endl;
				Eigen::VectorXf op2 = data_.col(j);
				Eigen::MatrixXf op3 = op1*op2.transpose();
				std::cout << op3<< std::endl;
				for (int k = 0; k < data_.rows(); ++k)
				{

				}
			}
		}
	}
	std::cout << "map x: " << map_x << ", map y: " << map_y << std::endl;
}

Som::~Som()
{
}

void Som::DetermineRadiusInitial()
{
	int ms = std::max(map_x, map_y);
	float map_radius_ini = fmaxf(1.0f, ((float)ms / 4));
	float map_radius_fin = fmaxf(1, (map_radius_ini / 4));
	int train_len = std::max(1, (int)(10 * ((map_x * map_y) / data_.rows())));
}

void Som::NInv(int n, int &height, int &width)
{
	if (n > codebook_->GetWeights().rows())
		return;
	width = (n - 1) / map_x + 1;
	height = n - (width - 1) * map_x;
}

void Som::CalculateUMatrix()
{
	int y = map_x;
	int x = map_y;
	int ux = 2 * x - 1;
	int uy = 2 * y - 1;

	Eigen::MatrixXf umat = Eigen::MatrixXf::Zero(uy, ux);
	typedef Eigen::Matrix<Eigen::VectorXf, Eigen::Dynamic, Eigen::Dynamic> Neurons;
	Neurons neurons(map_x, map_y);
	//Eigen::MatrixXf neurons(map_x, map_y);
	
	int count = 0;
	std::cout << "neurons" << std::endl;
	for (int j = 0; j < map_y; ++j)
	{
		for (int i = 0; i < map_x; ++i)
		{
			neurons(i, j) = codebook_->GetWeights().row(count);
			count++;
		}
	}
	IO *io = new IO();
	Eigen::MatrixXf codebook = codebook_->GetWeights();
	io->SaveMatrix(codebook, "codebook");

	for (int j = 0; j < y; ++j)
	{
		for (int i = 0; i < x; ++i)
		{
			if (i < (x - 1)) {
				float d1 = algorithm_->CalculateNeuronDistance(neurons(j, i), neurons(j, i+1));
				umat(2*j,2*i+1) = d1;
			}
			if (j < (y - 1)) {
				float d1 = algorithm_->CalculateNeuronDistance(neurons(j, i), neurons(j+1, i));
				umat(2*j+1, 2*i) = d1;
			}
			if (i < (x - 1) && j < (y - 1))
			{
				float d1 = algorithm_->CalculateNeuronDistance(neurons(j, i), neurons(j + 1, i + 1));
				float d2 = algorithm_->CalculateNeuronDistance(neurons(j + 1, i), neurons(j, i + 1));
				umat(2 * j + 1, 2 * i + 1) = (d1 + d2) / (2 * sqrt(2));
			}
		}
	}

	for (int j = 0; j < uy; j += 2)
	{
		for (int i = 0; i < ux; i += 2)
		{
			std::vector<float> a;
			if (i > 0 && j > 0 && i < (ux - 1) && j < (uy - 1)) // middle part of the map
			{
				a.push_back(umat(j, i - 1));
				a.push_back(umat(j, i + 1));
				a.push_back(umat(j - 1, i));
				a.push_back(umat(j + 1, i));
			}
			else if (j == 0 && i > 0 && i < (ux - 1)) // upper edge
			{
				a.push_back(umat(j, i - 1));
				a.push_back(umat(j, i + 1));
				a.push_back(umat(j + 1, i));
			}
			else if (j == (uy - 1) && i > 0 && i < (ux - 1)) // lower edge
			{
				a.push_back(umat(j, i - 1));
				a.push_back(umat(j, i + 1));
				a.push_back(umat(j - 1, i));
			}
			else if (i == 0 && j > 0 && j < (uy - 1)) // left edge
			{
				a.push_back(umat(j, i + 1));
				a.push_back(umat(j - 1, i));
				a.push_back(umat(j + 1, i));
			}
			else if (i == (ux - 1) && j > 0 && j < (uy - 1)) // right edge
			{
				a.push_back(umat(j, i - 1));
				a.push_back(umat(j - 1, i));
				a.push_back(umat(j + 1, i));
			}
			else if (i == 0 && j == 0) // top left corner
			{
				a.push_back(umat(j, i + 1));
				a.push_back(umat(j + 1, i));
			}
			else if (i == (ux - 1) && j == 0) // top right corner
			{
				a.push_back(umat(j, i - 1));
				a.push_back(umat(j + 1, i));
			}
			else if (i == 0 && j == (uy - 1)) // bottom left corner
			{
				a.push_back(umat(j, i + 1));
				a.push_back(umat(j - 1, i));
			}
			else if (i == (ux - 1) && j == (uy - 1)) // bottom right corner
			{
				a.push_back(umat(j, i - 1));
				a.push_back(umat(j - 1, i));
			}
			else
			{
				a.push_back(0.0f);
			}
			umat(j, i) = algorithm_->CalculateMedian(a);
		}

	}
	std::cout << "UMAT:" << std::endl << umat << std::endl;
	umat_ = umat;
	io->SaveUMAT(umat);
	delete io;
	/*Watershed * w = new Watershed();
	Eigen::MatrixXf r = w->transform(umat_);
	std::cout << "r" << std::endl << r;*/
}

void Som::CalculatePMatrix()
{
	ParetoDensity pd;
	Eigen::MatrixXf weigths = codebook_->GetWeights();
	Eigen::VectorXf resultP = pd.CalculateDensity(data_, weigths, 0.2);
	Eigen::MatrixXf PMatrix = algorithm_->Reshape(resultP, map_x, map_y);
	IO *io = new IO();
	io->SaveMatrix(PMatrix, "pmatrix");
}
