#include "Train.h"

#include <Eigen/Core>
#include <iostream>
#include <omp.h>

Train::Train(Som &som, Train::TrainType train_type, bool fine_adjustment) : fine_adjustment_(fine_adjustment)
{
	switch (train_type)
	{
	case TrainType::CLASSIC:
		ClassicSomTrain(som);
	default:
		break;
	}
}

bool Train::ClassicSomTrain(Som & som)
{
	std::cout << "Training... " << std::endl;
	float sigma = som.params_->sigma_;
	int dim = som.data_.cols();
	float dist, learning_rate = 0.0f;
	float hci_exp = 0.0f;
	Eigen::MatrixXf::Index index;
	Eigen::MatrixXf weights = som.codebook_->GetWeights();
	#pragma omp parallel for
	for (int current_epoch = 0; current_epoch < som.params_->train_len_; ++current_epoch)
	{
		sigma = som.algorithm_->Radius(som.params_->sigma_, current_epoch, som.params_->time_constant_);
		learning_rate = som.algorithm_->LearningRate(som.params_->learning_rate_, current_epoch);
		int rand_sample = rand() % (som.data_.rows() + 1);
		auto sample = som.data_.row(rand_sample);
		//get BMU
		(weights.rowwise() - sample).rowwise().squaredNorm().minCoeff(&index);
		auto bmu = som.grid_.row(index);
		for (int n = 0; n < weights.rows(); ++n)
		{
			dist = (bmu - som.grid_.row(n)).squaredNorm();
			float sigma2 = sigma*sigma;
			if (dist < sigma2)
			{
				hci_exp = std::exp((-dist) / (2 * sigma2));
				weights.row(n) += learning_rate * hci_exp * (sample - weights.row(n));
			}

		}
	}
	if (fine_adjustment_)
	{
		//Ajuste fino
		som.params_->SetSigma(10.0f);
		som.params_->SetLearningRate(0.1f);

		for (int current_epoch = 0; current_epoch < som.params_->train_len_; ++current_epoch)
		{
			sigma = som.algorithm_->Radius(som.params_->sigma_, current_epoch, som.params_->time_constant_);
			learning_rate = som.algorithm_->LearningRate(som.params_->learning_rate_, current_epoch);
			int rand_sample = rand() % (som.data_.rows() + 1);
			auto sample = som.data_.row(rand_sample);
			//get BMU
			(weights.rowwise() - sample).rowwise().squaredNorm().minCoeff(&index);
			auto bmu = som.grid_.row(index);
			for (int n = 0; n < weights.rows(); ++n)
			{
				dist = (bmu - som.grid_.row(n)).squaredNorm();
				float sigma2 = sigma*sigma;
				if (dist < sigma2)
				{
					hci_exp = std::exp((-dist) / (2 * sigma2));
					weights.row(n) += learning_rate * hci_exp * (sample - weights.row(n));
				}

			}
		}
	}

	som.codebook_->SetWeightsEndTrain(weights);
	std::cout << "Final weights" << std::endl;
	std::cout << weights << std::endl;

	return true;
}

bool Train::ClassicSomTrain2(Som & som)
{
	std::cout << "Training... " << std::endl;
	float sigma = som.params_->sigma_;
	int dim = som.data_.cols();
	float dist, learning_rate = 0.0f;
	float hci_exp = 0.0f;
	Eigen::MatrixXf::Index index;
	Eigen::MatrixXf weights = som.codebook_->GetWeights();
	#pragma omp parallel for
	for (int current_epoch = 0; current_epoch < som.params_->train_len_; ++current_epoch)
	{
		sigma = som.algorithm_->Radius(som.params_->sigma_, current_epoch, som.params_->time_constant_);
		learning_rate = som.algorithm_->LearningRate(som.params_->learning_rate_, current_epoch);
		for (int d = 0; d < som.data_.rows(); ++d) {
			(weights.rowwise() - som.data_.row(d)).rowwise().squaredNorm().minCoeff(&index);
			auto bmu = som.grid_.row(index);
			for (int n = 0; n < weights.rows(); ++n)
			{
				dist = (bmu - som.grid_.row(n)).squaredNorm();
				float sigma2 = sigma*sigma;
				//std::cout << n << std::endl;
				if (dist < sigma2)
				{
					hci_exp = std::exp((-dist) / (2 * sigma2));
					weights.row(n) += learning_rate * hci_exp * (som.data_.row(d) - weights.row(n));
				}
			}
		}
		std::cout << current_epoch << std::endl;
	}
	if (fine_adjustment_)
	{
		//Ajuste fino
		som.params_->SetSigma(10.0f);
		som.params_->SetLearningRate(0.1f);

		for (int current_epoch = 0; current_epoch < som.params_->train_len_; ++current_epoch)
		{
			sigma = som.algorithm_->Radius(som.params_->sigma_, current_epoch, som.params_->time_constant_);
			learning_rate = som.algorithm_->LearningRate(som.params_->learning_rate_, current_epoch);
			int rand_sample = rand() % (som.data_.rows() + 1);
			auto sample = som.data_.row(rand_sample);
			//get BMU
			(weights.rowwise() - sample).rowwise().squaredNorm().minCoeff(&index);
			auto bmu = som.grid_.row(index);
			for (int n = 0; n < weights.rows(); ++n)
			{
				dist = (bmu - som.grid_.row(n)).squaredNorm();
				float sigma2 = sigma*sigma;
				if (dist < sigma2)
				{
					hci_exp = std::exp((-dist) / (2 * sigma2));
					weights.row(n) += learning_rate * hci_exp * (sample - weights.row(n));
				}

			}
		}
	}

	som.codebook_->SetWeightsEndTrain(weights);
	std::cout << "Final weights" << std::endl;
	std::cout << weights << std::endl;

	return true;
}

Train::~Train()
{
}
