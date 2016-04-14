#include "Som.h"
#include <iostream>
#include <limits>
#include <vector>

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
	Eigen::MatrixXf umat(2 * map_x - 1, 2 * map_y - 1);
	typedef Eigen::Matrix<Eigen::VectorXf, Eigen::Dynamic, Eigen::Dynamic> Neurons;
	Neurons neurons(map_x, map_y);
	
	int count = 0;
	for (int i = 0; i < map_x; ++i)
	{
		for (int j = 0; j < map_y; ++j)
		{
			neurons(i, j) = codebook_->GetWeights().row(count);
			count++;
		}
	}

	for (int i = 0; i < map_x; ++i)
	{
		for (int j = 0; j < map_y; ++j)
		{
			if (j < (map_y - 1))
				umat(2*i,2*j+1) = algorithm_->CalculateNeuronDistance(neurons(i, j), neurons(i, j+1));
			if (i < (map_x - 1))
				umat(2*i+1, 2*j) = algorithm_->CalculateNeuronDistance(neurons(i, j), neurons(i+1, j));
			if (i < (map_x - 1) && j < (map_y - 1))
			{
				float d1 = algorithm_->CalculateNeuronDistance(neurons(i, j), neurons(i + 1, j + 1));
				float d2 = algorithm_->CalculateNeuronDistance(neurons(i + 1, j), neurons(i, j + 1));
				umat(2 * i + 1, 2 * j + 1) = (d1 + d2) / (2 * sqrt(2));
			}
		}
		std::vector<float> a;
		for (int i = 0; i < umat.rows(); i += 2)
		{
			for (int j = 0; j < umat.cols(); j += 2)
			{
				if (j > 0 && i > 0 && i < (umat.rows() - 1) && j < (umat.cols() - 1))
				{
					a.push_back(umat(i, j - 1));
					a.push_back(umat(i, j + 1));
					a.push_back(umat(i - 1, j));
					a.push_back(umat(i + 1, j));
				}
				else if (i == 0 && j > 0 && j < (umat.cols() - 1))
				{
					a.push_back(umat(i, j - 1));
					a.push_back(umat(i, j + 1));
					a.push_back(umat(i + 1, j));
				}
				else if (i == (umat.rows() - 1) && j > 0 && j < (umat.cols() - 1))
				{
					a.push_back(umat(i, j - 1));
					a.push_back(umat(i, j + 1));
					a.push_back(umat(i - 1, j));
				}
				else if (j == 0 && i > 0 && i < (umat.rows() - 1))
				{
					a.push_back(umat(i, j + 1));
					a.push_back(umat(i - 1, j));
					a.push_back(umat(i + 1, j));
				}
				else if (j == (umat.cols() - 1) && i > 0 && i < (umat.rows() - 1))
				{
					a.push_back(umat(i, j - 1));
					a.push_back(umat(i - 1, j));
					a.push_back(umat(i + 1, j));
				}
				else if (j == 0 && i == 0)
				{
					a.push_back(umat(i, j + 1));
					a.push_back(umat(i + 1, j));
				}
				else if (i == 0 && j > 0 && j < (umat.cols() - 1))
				{
					a.push_back(umat(i, j));
					a.push_back(umat(i, j));
					a.push_back(umat(i, j));
				}
				else if (j == (umat.rows() - 1) && i == 0)
				{
					a.push_back(umat(i, j - 1));
					a.push_back(umat(i + 1, j));
				}
				else if (j == 0 && i < (umat.rows() - 1))
				{
					a.push_back(umat(i, j + 1));
					a.push_back(umat(i - 1, j));
				}
				else if (j == (umat.cols() - 1) && i == (umat.rows() - 1))
				{
					a.push_back(umat(i, j - 1));
					a.push_back(umat(i - 1, j));
				}
				else
				{
					a.push_back(0.0f);
				}
				umat(i, j) = algorithm_->CalculateMedian(a);
				//a.clear();
			}

		}
	}
	std::cout << "UMAT:" << std::endl << umat << std::endl;
}
