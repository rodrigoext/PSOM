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
	InitGrid(topology);
	Train();
	CalculateUMatrix();
	CalculatePMatrix();
}

void Som::InitGrid(Som::Topology topology)
{
	grid_ = Eigen::MatrixXf(codebook_->GetWeights().rows(), 2);
	int width, height;
	switch (topology)
	{
	case HEXAGONAL:
		for (int i = 0; i <= map_x*map_y; ++i)
		{
			NInv(i, height, width);
			grid_(i, 0) = static_cast<float>(width);
			grid_(i, 1) = static_cast<float>(height) * 0.8660;
			if (height % 2 == 0)
				grid_(i, 0) += 0.5;
		}
		break;
	case RETANGULAR:
		for (int i = 0; i < map_x*map_y; ++i)
		{
			/*NInv(i, height, width);
			grid_(i - 1, 0) = static_cast<float>(width);
			grid_(i - 1, 1) = static_cast<float>(height);*/
			NInv(i, height, width);
			grid_(i, 0) = static_cast<float>(width);
			grid_(i, 1) = static_cast<float>(height);
		}
		break;
	default:
		break;
	}
}

void Som::Train()
{
	std::cout << "Training... ";
	float sigma = params_->sigma_;
	int dim = data_.cols();
	float dist, learning_rate = 0.0f;
	float hci_exp = 0.0f;
	Eigen::MatrixXf::Index index;
	Eigen::MatrixXf weights = codebook_->GetWeights();
	//#pragma omp parallel for
	for (int current_epoch = 0; current_epoch < params_->train_len_; ++current_epoch)
	{
		sigma = algorithm_->Radius(params_->sigma_, current_epoch, params_->time_constant_);
		learning_rate = algorithm_->LearningRate(params_->learning_rate_, current_epoch);
		int rand_sample = rand() % (data_.rows());
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

	if (params_->fine_tune_)
	{
		std::cout << "ajuste fino" << std::endl;
		//Fine tune
		sigma = 0.1f;
		learning_rate = 0.01f;

		for (int current_epoch = 0; current_epoch < 5*params_->train_len_; ++current_epoch)
		{
			int rand_sample = rand() % (data_.rows());
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
	}
	codebook_->SetWeightsEndTrain(weights);
	std::cout << "OK!" << std::endl;
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
	if (n >= codebook_->GetWeights().rows())
		return;
	width = n / map_x;
	height = n - width * map_x;
}

void Som::CalculateUMatrix()
{
	IO *io = new IO();
	//Load for test
	//Eigen::MatrixXf codebook_load = io->LoadCSV("codebook_bom");
	//codebook_->SetWeightsEndTrain(codebook_load);
	int y = map_x;
	int x = map_y;
	int ux = 2 * x - 1;
	int uy = 2 * y - 1;

	Eigen::MatrixXf umat = Eigen::MatrixXf::Zero(uy, ux);
	typedef Eigen::Matrix<Eigen::VectorXf, Eigen::Dynamic, Eigen::Dynamic> Neurons;
	Neurons neurons(map_x, map_y);
	//Eigen::MatrixXf neurons(map_x, map_y);

	int count = 0;
	for (int j = 0; j < map_y; ++j)
	{
		for (int i = 0; i < map_x; ++i)
		{
			neurons(i, j) = codebook_->GetWeights().row(count);
			count++;
		}
	}
	
	Eigen::MatrixXf codebook = codebook_->GetWeights();
	io->SaveMatrix(codebook, "codebook");
	std::cout << "Calculating U-Matrix... ";
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
	std::cout << "OK!" << std::endl;
	umat_ = algorithm_->FilterMedian(umat);
	io->SaveUMAT(umat);
	io->SaveMatrix(umat_, "umatfiltred");
	delete io;
}

void Som::CalculatePMatrix()
{
	std::cout << "Calculating PMatrix" << std::endl;
	ParetoDensity pd;
	Eigen::MatrixXf weigths = codebook_->GetWeights();
	Eigen::VectorXf resultP = pd.CalculateDensity(data_, weigths, umat_.maxCoeff());
	Eigen::MatrixXf p_matrix = algorithm_->Reshape(resultP, map_x, map_y);
	//pmat_ = algorithm_->FilterMedian(p_matrix);
	pmat_ = p_matrix;
	IO *io = new IO();
	io->SaveMatrix(pmat_, "pmatrix");
	io->SaveMatrix(data_, "data");
	delete io;
	std::cout << pmat_ << std::endl;
}

void Som::CalculateAllMatrix() {
	std::cout << "Calculating..." << std::endl;
	IO * io = new IO();
	Algorithm * algo = new Algorithm();

	//Load for test
    //Eigen::MatrixXf ustar_load = io->LoadCSV("ustar_bom");

	Eigen::MatrixXf umu = CalculateUMatrixUltsch();
	Eigen::MatrixXf ustar = CalculateUStarMatrix(umu, pmat_) * 10;
	io->SaveMatrix(ustar, "ustar");
	Watershed *w = new Watershed();
    //Eigen::MatrixXf ustar_w2 = w->transform_v2(algo->FilterMedian(ustar));
    Eigen::MatrixXf ustar_w = w->transform_v2(algo->FilterMedian(ustar));
	//ustar_w = w->transform(ustar_load);
	//Eigen::MatrixXf um_w = w->transform(umu);
	io->SaveMatrix(ustar_w, "ustarw");
	//io->SaveMatrix(um_w, "umatw");
	//ustar_w = io->LoadCSV("ustarw_bom");
	std::cout << "Calculating Immersion" << std::endl;
	Eigen::MatrixXf imm = CalculateImmersion(pmat_, ustar_w);
	io->SaveMatrix(imm, "immersion");
	std::cout << "Simulating" << std::endl;
	Eigen::VectorXf sim = SimulateClustering(data_, ustar_w, imm);
	io->SaveVector(sim, "simulation");
	Eigen::VectorXf simP = SimulateClusteringParallel(data_, ustar_w, imm);
	io->SaveVector(simP, "simulationP");
	delete io;
	delete w;
}

Eigen::MatrixXf Som::CalculateUStarMatrix(Eigen::MatrixXf &umat, Eigen::MatrixXf &pmat)
{
	std::cout << "Calculating U*-Matrix... ";
	int linhas = umat.rows();
	int colunas = umat.cols();

	std::cout << linhas << " " << colunas << std::endl;
	Eigen::MatrixXf ustarmat(linhas,colunas);

	float plow;

	float meanp,maxp,maxp2;
	meanp = pmat.sum() / (linhas*colunas);
	maxp = pmat.maxCoeff();
	for (int l = 0; l < linhas; l++){
		for (int c = 0; c < colunas; c++){
			plow = CalculatePlow(pmat,l,c);
			ustarmat(l,c) = umat(l,c) * plow;
		}
	}
	std::cout << "OK!" << std::endl;
	IO * io = new IO();
	io->SaveMatrix(ustarmat, "ustar");
	delete io;
	return ustarmat;
}

float Som::CalculatePlow(Eigen::MatrixXf &pmat, int li, int ci)
{
	int linhas = pmat.rows();
	int colunas = pmat.cols();
	int counter = 0;
	for (int l = 0; l < linhas; l++) {
		for (int c = 0; c < colunas; c++) {
			if (pmat(li,ci) < pmat(l,c))
				counter++;
		}
	}
	return (float)counter / (float)(linhas*colunas);
}

Eigen::MatrixXf Som::CalculateUMatrixUltsch()
{
	std::cout << "Calculating U-Matrix Ultsch... ";
	Eigen::MatrixXf umatriz(map_x, map_y);
	int na, nl,counter;
	float du,di;
	int linha,coluna;
	float max = 0;
	float dtemp;
	int unidades = map_x*map_y;
	for (int i = 0; i < unidades; i++) {
		NInv(i, na, nl);
		counter = 0;
		du = 0;

		for (int v1 = -1; v1 < 2; v1++){
			if (v1==0)
				v1++;
			linha = na + v1;
			if (linha >= 0 and linha < map_x){
				//dtemp = d(i,n(linha,nl));
				dtemp = (codebook_->GetWeights().row(i) - codebook_->GetWeights().row(linha+nl*map_x)).squaredNorm();
				du += dtemp;

				if (max < dtemp)
					max = dtemp;
				counter++ ;
			}
		}

		for (int v1 = -1; v1 < 2; v1++){
			if (v1==0)
				v1++;
			coluna = nl + v1;
			if (coluna >= 0 and coluna < map_y){
				//dtemp = d(i,n(na,coluna));
				dtemp = (codebook_->GetWeights().row(i) - codebook_->GetWeights().row(na+coluna*map_x)).squaredNorm();
				du += dtemp;

				if (max < dtemp)
					max = dtemp;
				counter++ ;
			}
		}
		umatriz(na,nl) = (du / counter);
	}
	std::cout << "OK!" << std::endl;
	IO * io = new IO();
	io->SaveMatrix(umatriz, "umultsch");
	delete io;
	return umatriz;
}

Eigen::MatrixXf Som::CalculateImmersion(Eigen::MatrixXf &pmat, Eigen::MatrixXf &umat)
{
	Eigen::MatrixXf result(pmat.rows(),pmat.cols());
	Eigen::MatrixXf umattemp = umat*-1;
	int temp = 0;
	int ltemp,ctemp;
	for (int l = 0 ; l < map_x ; l++){
		for (int c = 0; c < map_y; c++){
			temp = CalculateImersion(l,c,umattemp);
			NInv(temp,ltemp,ctemp);
			result(l,c) = CalculateImersion(ltemp,ctemp,pmat);
		}
	}
	return result;
}

int Som::CalculateImersion(int linha, int coluna, Eigen::MatrixXf &mat) {
	double max;

	int maxIt = (map_x*map_y)/2;
	int itAtual = 0;

	int linhaAtual = linha,colunaAtual = coluna;
	int linhaFinal = 0,colunaFinal = 0;
	int lverif,cverif;

	max = mat(linha,coluna);

	while (true) {
		lverif = linhaAtual - 1;
		if(lverif >= 0){
			if(mat(linhaAtual,colunaAtual) <= mat(lverif,colunaAtual)){
				max = mat(lverif,colunaAtual);
				linhaFinal = lverif;
				colunaFinal = colunaAtual;
			}
		}

		lverif = linhaAtual + 1;
		if(lverif < map_x){
			if(mat(linhaAtual,colunaAtual) <= mat(lverif,colunaAtual))
				if (mat(lverif,colunaAtual) >= max){
					max = mat(lverif,colunaAtual);
					linhaFinal = lverif;
					colunaFinal = colunaAtual;
				}
		}

		cverif = colunaAtual - 1;
		if(cverif >= 0){
			if(mat(linhaAtual,colunaAtual) <= mat(linhaAtual,cverif))
				if (mat(linhaAtual,cverif) >= max){
					max = mat(linhaAtual,cverif);
					linhaFinal = linhaAtual;
					colunaFinal = cverif;
				}
		}

		cverif = colunaAtual + 1;
		if(cverif < map_y){
			if(mat(linhaAtual,colunaAtual) <= mat(linhaAtual,cverif))
				if (mat(linhaAtual,cverif) >= max){
					max = mat(linhaAtual,cverif);
					linhaFinal = linhaAtual;
					colunaFinal = cverif;
				}
		}

		if (linhaAtual==linhaFinal and colunaAtual==colunaFinal)
			break;

		if (linhaFinal==0 or colunaFinal==0){
			linhaFinal = linhaAtual;
			colunaFinal = colunaAtual;
			break;
		}

		linhaAtual = linhaFinal;
		colunaAtual = colunaFinal;

		if (itAtual++ > maxIt)
			break;
	}

	return linhaFinal+colunaFinal*map_x;
}

Eigen::VectorXf Som::SimulateClustering(Eigen::MatrixXf &data, Eigen::MatrixXf &watershed, Eigen::MatrixXf &immersion)
{
	Eigen::VectorXf result(data.rows());

	float min, dist;
	int menorN;
	int lin,col;
	int tempMap;
	int temp1, temp2;
	temp1 = 1;
//#pragma omp parallel for
	for (int i = 0 ; i < data.rows() ; i++){
		min = algorithm_->CalculateNeuronDistance(codebook_->GetWeights().row(1), data.row(i));
		menorN = 0;
		for (int j = 0; j < map_x*map_y ; j++ ){
			dist = algorithm_->CalculateNeuronDistance(codebook_->GetWeights().row(j), data.row(i));
			if (dist < min){
				menorN = j;
				min = dist;
			}
		}
		NInv(menorN,lin,col);
		tempMap = immersion(lin,col);
		NInv(tempMap,lin,col);
		temp1 = watershed(lin,col);
		if(temp1 != -2)
			temp2 = temp1;
		result(i) = temp2;
	}

	return result;
}
Eigen::VectorXf Som::SimulateClusteringParallel(Eigen::MatrixXf &data, Eigen::MatrixXf &watershed, Eigen::MatrixXf &immersion)
{
	Eigen::VectorXf result(data.rows());

	float min, dist;
	int menorN;
	int lin,col;
	int tempMap;
	int temp1, temp2;
	temp1 = 1;
	#pragma omp parallel shared(result)
	for (int i = 0 ; i < data.rows() ; i++){
		min = algorithm_->CalculateNeuronDistance(codebook_->GetWeights().row(1), data.row(i));
		menorN = 0;
		for (int j = 0; j < map_x*map_y ; j++ ){
			dist = algorithm_->CalculateNeuronDistance(codebook_->GetWeights().row(j), data.row(i));
			if (dist < min){
				menorN = j;
				min = dist;
			}
		}
		NInv(menorN,lin,col);
		tempMap = immersion(lin,col);
		NInv(tempMap,lin,col);
		temp1 = watershed(lin,col);
		if(temp1 != -2)
			temp2 = temp1;
		result(i) = temp2;
	}

	return result;
}
