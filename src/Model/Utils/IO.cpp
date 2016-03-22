#include "IO.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace Eigen;

IO::IO()
{
}

MatrixXf IO::LoadData(const char * file_name)
{
	vector<vector<float>> data_set = ReadData(file_name);
	
	MatrixXf data_e = Vector2EingenMatrix(data_set);

	return data_e;
}

vector<vector<float>> IO::ReadData(const char * file_name)
{
	vector<vector<string>> data;
	vector<vector<float>> data_float;

	ifstream infile;

	infile.open(file_name);

	cerr << "Reading data..." << endl;

	int i, j = 0;

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

MatrixXf IO::Vector2EingenMatrix(vector<vector<float>> d)
{
	MatrixXf data_eigen(d.size(), d[0].size());

	for (int i = 0; i < d.size(); ++i)
	{
		for (int j = 0; j < d[0].size(); ++j)
		{
			data_eigen(i, j) = d[i][j];
		}
	}

	return data_eigen;
}

IO::~IO()
{
}
