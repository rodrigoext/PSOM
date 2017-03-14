#pragma once

#include "NeuralNetwork.h"
#include "Node.h"
#include "Model/Utils/Parameter.h"
#include "Model/Utils/Algorithm.h"
#include "Model/Network/Codebook.h"

#include <Eigen>
#include <memory>
#include <map>
class Som : public NeuralNetwork
{
protected:
	Eigen::MatrixXf data_;
	Eigen::MatrixXf grid_;
	Eigen::MatrixXf umat_;
	Eigen::MatrixXf pmat_;
    Eigen::MatrixXf ustarmat_;
    Eigen::MatrixXf imm_;
    Eigen::MatrixXf ustarw_;
    Eigen::VectorXi bmu_;
    Eigen::VectorXf class_;
	std::shared_ptr<Parameter> params_;
	std::shared_ptr<Codebook> codebook_;
	std::shared_ptr<Algorithm> algorithm_;
    int map_x, map_y;
    bool aaaaa_;

public:
	enum Topology
	{
		RETANGULAR,
		HEXAGONAL
	};
	Som(Eigen::MatrixXf data, Som::Topology topology = Topology::RETANGULAR);
	Som(Eigen::MatrixXf data, std::shared_ptr<Parameter> params, Som::Topology topology = Topology::RETANGULAR);
	void Train();
	Eigen::MatrixXf GetUMatrix() {return umat_;}
    Eigen::MatrixXf GetUStarMatrix() {return ustarmat_;}
    Eigen::MatrixXf GetUStarMatrixWatershed() {return ustarw_;}
	Eigen::MatrixXf CalculateImmersion(Eigen::MatrixXf &pmat, Eigen::MatrixXf &umat);
	Eigen::VectorXf SimulateClustering(Eigen::MatrixXf &data, Eigen::MatrixXf &watershed, Eigen::MatrixXf &immersion);
	Eigen::VectorXf SimulateClusteringParallel(Eigen::MatrixXf &data, Eigen::MatrixXf &watershed, Eigen::MatrixXf &immersion);
    Eigen::VectorXf SimulateWithNeurons();
    void Simulate();
    void ClusterMap();
	void CalculateAllMatrix();
    void UpdateDiscretizarion(int value);
    Eigen::VectorXf GetClassification() { return class_; }
	virtual ~Som();
private:
	void InitGrid(Som::Topology topology);
	void CalculateMapSize();
	void DetermineRadiusInitial();
	void NInv(int n, int &width, int &height);
	void CalculateUMatrix();
	void CalculatePMatrix();
	Eigen::MatrixXf CalculateUStarMatrix(Eigen::MatrixXf &umat, Eigen::MatrixXf &pmat);
	float CalculatePlow(Eigen::MatrixXf &pmat, int li, int ci);
	Eigen::MatrixXf CalculateUMatrixUltsch();
	int CalculateImersion(int linha, int coluna, Eigen::MatrixXf &mat);
	//Node n;
};
