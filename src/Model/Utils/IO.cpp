#include "IO.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Eigen;

std::string path_to_save = "/home/rodrigo/l3c/plugins-v3o2/PluginL3C/libs/ltrace/kohonen/src/Data/";

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");

IO::IO()
{
}

MatrixXf IO::LoadData(const char * file_name, bool normalize)
{
	vector<vector<float> > data_set = ReadData(file_name);

	MatrixXf data_e = Vector2EingenMatrix(data_set);
	if (normalize)
		MapMinMax(data_e);
	std::cout << " OK!" << std::endl;
	return data_e;
}

vector<vector<float> > IO::ReadData(const char * file_name)
{
	vector<vector<string> > data;
	vector<vector<float> > data_float;

	ifstream infile;

	infile.open(file_name);

	cerr << "Reading data...";

	if (infile.is_open())
	{
		while (infile)
		{
			string s;

			if (!getline(infile, s))
				break;

			istringstream ss(s);
			vector<float> record;

			while (ss)
			{
				string s2;
				if (!getline(ss, s2, ','))
					break;

				float temp = 0;

				istringstream iss(s2);

				iss >> temp;

				record.push_back(temp);

			}

			data_float.push_back(record);
		}
	}

	if (!infile.eof())
	{
		cerr << "Read data -> OK!\n";
	}

	return data_float;
}

void IO::SaveUMAT(Eigen::MatrixXf &data)
{
    std::ofstream file(path_to_save+"umat.csv");
	if (file.is_open())
	{
		file << data.format(CSVFormat);
	}
}

void IO::SaveMatrix(Eigen::MatrixXf &data, std::string file_name) {
    std::ofstream file(path_to_save+file_name+".csv");
	if (file.is_open())
	{
		file << data.format(CSVFormat);
	}
}

void IO::SaveVector(Eigen::VectorXf &data, std::string file_name)
{
    std::ofstream file(path_to_save+file_name+".csv");
	if (file.is_open())
	{
		file << data.format(CSVFormat);
	}
}

void IO::SaveVectorINT(Eigen::VectorXi &data, std::string file_name)
{
    std::ofstream file(path_to_save+file_name+".csv");
    if (file.is_open())
    {
        file << data.format(CSVFormat);
    }
}

void IO::MapMinMax(Eigen::MatrixXf &data)
{
	float min = 0.0f;
	float max = 0.0f;

	float scale_range = 1.0f - (-1.0f);
	

	for (int i = 0; i < data.cols(); ++i)
	{
		min = data.col(i).minCoeff();
		max = data.col(i).maxCoeff();

		float value_range = max - min;

		for (int j = 0; j < data.col(i).size(); ++j)
		{
			data(j, i) = (data(j, i) - min) / value_range * scale_range + (-1.0f);
			//data[i][j] = (data[i][j] - min) / value_range * scale_range + (-1.0f);
		}
	}
}

MatrixXf IO::Vector2EingenMatrix(vector<vector<float> > d)
{
	MatrixXf data_eigen(d.size(), d[0].size());

	for (unsigned int i = 0; i < d.size(); ++i)
	{
		for (unsigned int j = 0; j < d[0].size(); ++j)
		{
			data_eigen(i, j) = d[i][j];
		}
	}

	return data_eigen;
}

Eigen::MatrixXf IO::LoadCSV(std::string file_name) {
	std::ifstream indata;
    indata.open(path_to_save+file_name+".csv");
    std::string line;
    std::vector<float> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stof(cell));
        }
        ++rows;
    }
	Eigen::MatrixXf retorno = Eigen::Map<Eigen::MatrixXf>(values.data(), rows, values.size()/rows); 
    return retorno;
}

IO::~IO()
{
}
