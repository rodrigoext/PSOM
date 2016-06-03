/*
 * Watershed.cpp
 *
 *  Created on: 15/05/2016
 *      Author: rodrigo
 */

#include <Model/Utils/Watershed.h>
#include <queue>
#include <iostream>
// Níveis de discretização da imagem
#define LEVELS 31
# define MASK -2 /* initial value of a threshold level */
# define WSHED 0 /* value of the pixels belonging to the watersheds */
# define INIT -1 /* initial value of im, */

Watershed::Watershed() {
	// TODO Auto-generated constructor stub

}

Watershed::~Watershed() {
	// TODO Auto-generated destructor stub
}

Eigen::MatrixXf Watershed::transform(const Eigen::MatrixXf &input) {
	Eigen::MatrixXf imi = input;
		int rows = imi.rows();
		int cols = imi.cols();
		int units = rows * cols;
		imi = imi * LEVELS;

		for (int l = 0; l < rows; l++) {
			for (int c = 0; c < cols; c++) {
				imi(l, c) = ((int) (imi(l, c) + 0.5)) + 1;
			}
		}

		std::queue<pos> fifo;

		Eigen::MatrixXf imo(imi);

		imo.setConstant(INIT);

		int current_label = 0;
		int current_distance;

		Eigen::MatrixXf imd(imo);
		imd.setConstant(0);

		Eigen::MatrixXf temp_mat(imi);
		Eigen::VectorXf sorted(rows*cols);

		std::vector<pos> positions;
		int row_min = 0, col_min = 0;
		float temp_min = imi(0,0);

		float hmax = temp_mat.maxCoeff();
		float hmin = temp_mat.minCoeff();

		for(int i = 0; i < units; ++i)
		{
			for (int r = 0; r < rows; ++r)
			{
				for (int c = 0; c < cols; ++c)
				{
					if (temp_min > temp_mat(r,c)) {
						temp_min = temp_mat(r,c);
						row_min = r;
						col_min = c;
					}
				}
			}
			pos p;
			p.row = row_min;
			p.col = col_min;
			positions.push_back(p);
			sorted(i) = temp_min;
			temp_mat(row_min, col_min) = hmax + 1;
			row_min = 0;
			col_min = 0;
			temp_min = temp_mat(0, 0);
		}
		int last_pixel = 0;
		pos position;
		for (int h = hmin; h < hmax; ++h)
		{
			while (h == sorted(last_pixel))
			{
				position = positions.at(last_pixel);
				imo(position.row, position.col) = MASK;
				std::vector<pos> positions_temp = vizinhanca4(position.row, position.col, rows, cols);
				for (unsigned int i = 0; i < positions_temp.size(); ++i)
				{
					pos p = positions_temp.at(i);
					if (imo(p.row, p.col) >= 0) {
						imd(p.row, p.col) = 1;
						fifo.push(p);
					}
				}
				last_pixel++;
				if (last_pixel > (units - 1))
					break;
			}
			current_distance = 1;
			pos *p = new pos;
			p->row = 0;
			fifo.push(*p);
			do {
				*p = fifo.front();
				fifo.pop();
				if (p->row == 0) {
					if (fifo.empty()) {
						break;
					} else {
						p = new pos;
						p->row = 0;
						fifo.push(*p);
						current_distance += 1;
						*p = fifo.front();
						fifo.pop(); // p = fifo.first();
					}
				}
				std::vector<pos> posicoes = vizinhanca4(p->row, p->col, rows, cols);
				for (unsigned int i = 0; i < posicoes.size(); i++) {
					pos pl = posicoes.at(i);

					// se plinha pertence a alguma bacia ou watershed
					if (imd(pl.row, pl.col) < current_distance and (imo(pl.row, pl.col) > 0
							or imo(pl.row, pl.col) == WSHED)) {

						if (imo(pl.row, pl.col) > 0) {
							if (imo(p->row, p->col) == MASK or imo(p->row, p->col) == WSHED)
								imo(p->row, p->col) = imo(pl.row, pl.col);
							else if (imo(p->row, p->col) != imo(pl.row, pl.col))
								imo(p->row, p->col) = WSHED;
						} else if (imo(p->row, p->col) == MASK)
							imo(p->row, p->col) = WSHED;
					} else if (imo(pl.row, pl.col) == MASK and imd(pl.row, pl.col) == 0) {
						imd(pl.row, pl.col) = current_distance + 1;
						fifo.push(pl);
					}
				}
			} while (1);

			int lastPixel2 = 0;
			double sortlast = sorted(lastPixel2);
			while (h == sorted(lastPixel2)) {
				position = positions.at(lastPixel2);
				imd(position.row, position.col) = 0;
				if (imo(position.row, position.col) == MASK) {
					current_label += 1;
					fifo.push(position);
					imo(position.row, position.col) = current_label;
					while (!fifo.empty()) {
						pos pl = fifo.front();
						fifo.pop(); // pl = fifo.first()
						std::vector<pos> posicoes = vizinhanca4(pl.row, pl.col, rows, cols);
						for (unsigned int i = 0; i < posicoes.size(); i++) {
							pos pll = posicoes.at(i);
							if (imo(pll.row, pll.col) == MASK) {
								fifo.push(pll);
								imo(pll.row, pll.col) = current_label;
							}
						}
					}
				}

				lastPixel2++;
				if (lastPixel2 > (units - 1))
					break;
			}
		}

		return imo;
}

// Funcao vizinhanca dos pixels com conectividade 4
std::vector<pos> Watershed::vizinhanca4(int lin, int col, int linhas, int colunas) {
	std::vector<pos> positions;
	if (col - 1 >= 0) {
		pos p;
		p.col = col - 1;
		p.row = lin;
		positions.push_back(p);
	}

	if (col + 1 < colunas) {
		pos p;
		p.col = col + 1;
		p.row = lin;
		positions.push_back(p);
	}
	if (lin - 1 >= 0) {
		pos p;
		p.col = col;
		p.row = lin - 1;
		positions.push_back(p);
	}
	if (lin + 1 < linhas) {
		pos p;
		p.col = col;
		p.row = lin + 1;
		positions.push_back(p);
	}
	return positions;
}

