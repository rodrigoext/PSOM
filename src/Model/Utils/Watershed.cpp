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

Watershed::Watershed() {
	// TODO Auto-generated constructor stub

}

Watershed::~Watershed() {
	// TODO Auto-generated destructor stub
}

Eigen::MatrixXf Watershed::transform(const Eigen::MatrixXf &input) {
	Eigen::MatrixXf imi = input;
	int linhas = imi.rows(), colunas = imi.cols();
	int unidades = linhas * colunas;
	for (int l = 0; l < linhas; l++) {
			for (int c = 0; c < colunas; c++) {
				imi(l, c) = imi(l, c) - imi.minCoeff();
			}
		}
	for (int l = 0; l < linhas; l++) {
			for (int c = 0; c < colunas; c++) {
				imi(l, c) = imi(l, c) / imi.maxCoeff();
			}
		}
	for (int l = 0; l < linhas; l++) {
			for (int c = 0; c < colunas; c++) {
				imi(l, c) = imi(l, c) *LEVELS;
			}
		}
	//imi = imi - imi.minCoeff();
	//imi = imi / imi.maxCoeff();
	//imi = imi * LEVELS;

	for (int l = 0; l < linhas; l++) {
		for (int c = 0; c < colunas; c++) {
			imi(l, c) = ((int) (imi(l, c) + 0.5)) + 1;
		}
	}

	std::queue<pos> fifo;

# define MASK -2 /* initial value of a threshold level */
# define WSHED 0 /* value of the pixels belonging to the watersheds */
# define INIT -1 /* initial value of im, */

	//. -input: imi, decimal image;
	Eigen::MatrixXf imo(imi); //-output: im0 image of the labeled watersheds;

	// Initializations:
	//-Value INIT is assigned to each pixel of im, :
	for (int i = 0; i < imo.rows(); ++i)
		for (int j = 0; j < imo.cols(); ++j)
			imo(i,j) = INIT;

	int currentLabel = 0;
	int currentDistance;
	Eigen::MatrixXf imd(imo);
	for (int i = 0; i < imo.rows(); ++i)
			for (int j = 0; j < imo.cols(); ++j)
				imd(i,j) = INIT;


	// Sort the pixels of imi in the increasing order of their gray values.
	Eigen::MatrixXf tempMat = imi;
	Eigen::VectorXf sorted(linhas * colunas);
	std::vector<pos> positions;

	int linhamin = 1, colunamin = 1;
	double tempMin = imi(1, 1);
	std::cout << "aaaqui" << std::endl;
	// Let hmin and hmax designate the lowest and highest values, respectively.
	double hmax = (double)tempMat.maxCoeff();
	double hmin = (double)tempMat.minCoeff();

	// Sort pixels
	for (int i = 0; i < linhas * colunas; i++) {
		for (int l = 0; l < linhas; l++) {
			for (int c = 0; c < colunas; c++) {
				if (tempMin > tempMat(l, c)) {
					tempMin = tempMat(l, c);
					linhamin = l;
					colunamin = c;
				}
			}
		}
		pos p;
		p.l = linhamin;
		p.c = colunamin;
		positions.push_back(p);
		sorted(i) = tempMin;
		tempMat(linhamin, colunamin) = hmax + 1;
		linhamin = 1;
		colunamin = 1;
		tempMin = tempMat(1, 1);
	}
	std::cout << "aaaqui" << std::endl;
	int lastPixel = 1;
	pos posicao;
	// For h = hmin to hmax {
	for (int h = hmin-1; h < hmax; h++) {
		// cout << "h: " << h << endl;
		// geodesic SKIZ of level h - 1 inside level h
		//For every pixel p such that imi(p) = h {
		while (h == sorted(lastPixel)) {
			// cout << "lastPixel: " << lastPixel << endl;
			posicao = positions.at(lastPixel - 1);
			std::cout << "aaaqui" << std::endl;
			imo(posicao.l, posicao.c) = MASK;
			std::vector<pos> posicoes = vizinhanca4(posicao.l, posicao.c, linhas,
					colunas);
			for (unsigned int i = 0; i < posicoes.size(); i++) {
				pos p = posicoes.at(i);
				if (imo(p.l, p.c) >= 0) {
					imd(p.l, p.c) = 1;
					fifo.push(p);
				}
			}
			lastPixel++;
			if (lastPixel > unidades)
				break;
		}
		currentDistance = 1;
		pos * p = new pos;
		p->l = 0;
		fifo.push(*p);
		do {
			*p = fifo.front();
			fifo.pop(); // p = fifo.first();
			if (p->l == 0) {
				if (fifo.empty())
					break;
				else {
					p = new pos;
					p->l = 0;
					fifo.push(*p);
					currentDistance += 1;
					*p = fifo.front();
					fifo.pop(); // p = fifo.first();
				}
			}
			std::vector<pos> posicoes = vizinhanca4(p->l, p->c, linhas, colunas);
			for (unsigned int i = 0; i < posicoes.size(); i++) {
				pos pl = posicoes.at(i);

				// se plinha pertence a alguma bacia ou watershed
				if (imd(pl.l, pl.c) < currentDistance and (imo(pl.l, pl.c) > 0
						or imo(pl.l, pl.c) == WSHED)) {

					if (imo(pl.l, pl.c) > 0) {
						if (imo(p->l, p->c) == MASK or imo(p->l, p->c) == WSHED)
							imo(p->l, p->c) = imo(pl.l, pl.c);
						else if (imo(p->l, p->c) != imo(pl.l, pl.c))
							imo(p->l, p->c) = WSHED;
					} else if (imo(p->l, p->c) == MASK)
						imo(p->l, p->c) = WSHED;
				} else if (imo(pl.l, pl.c) == MASK and imd(pl.l, pl.c) == 0) {
					imd(pl.l, pl.c) = currentDistance + 1;
					fifo.push(pl);
				}
			}

		} while (1);

		//For every pixel p such that imi(p) = h {
		int lastPixel2 = 1;
		double sortlast = sorted(lastPixel2);
		while (h == sorted(lastPixel2)) {
			posicao = positions.at(lastPixel2 - 1);
			imd(posicao.l, posicao.c) = 0;
			if (imo(posicao.l, posicao.c) == MASK) {
				currentLabel += 1;
				fifo.push(posicao);
				imo(posicao.l, posicao.c) = currentLabel;
				while (!fifo.empty()) {
					pos pl = fifo.front();
					fifo.pop(); // pl = fifo.first()
					std::vector<pos> posicoes = vizinhanca4(pl.l, pl.c, linhas,colunas);
					for (unsigned int i = 0; i < posicoes.size(); i++) {
						pos pll = posicoes.at(i);
						if (imo(pll.l, pll.c) == MASK) {
							fifo.push(pll);
							imo(pll.l, pll.c) = currentLabel;
						}
					}
				}
			}


			lastPixel2++;
			if (lastPixel2 > unidades)
				break;
		}
	}

	return imo;
}

// Funcao vizinhanca dos pixels com conectividade 4
std::vector<pos> Watershed::vizinhanca4(int lin, int col, int linhas, int colunas) {
	std::vector<pos> positions;
	if (col - 1 > 0) {
		pos p;
		p.c = col - 1;
		p.l = lin;
		positions.push_back(p);
	}

	if (col + 1 <= colunas) {
		pos p;
		p.c = col + 1;
		p.l = lin;
		positions.push_back(p);
	}
	if (lin - 1 > 0) {
		pos p;
		p.c = col;
		p.l = lin - 1;
		positions.push_back(p);
	}
	if (lin + 1 <= linhas) {
		pos p;
		p.c = col;
		p.l = lin + 1;
		positions.push_back(p);
	}
	/* conectividade8
	if (col - 1 > 0 and lin - 1 > 0) {
		pos p;
		p.c = col - 1;
		p.l = lin - 1;
		positions.push_back(p);
	}

	if (col + 1 <= colunas and lin + 1 <= linhas) {
		pos p;
		p.c = col + 1;
		p.l = lin + 1;
		positions.push_back(p);
	}
	if (lin - 1 > 0 and col + 1 <= colunas) {
		pos p;
		p.c = col + 1;
		p.l = lin - 1;
		positions.push_back(p);
	}
	if (lin + 1 <= linhas and col - 1 > 0) {
		pos p;
		p.c = col - 1;
		p.l = lin + 1;
		positions.push_back(p);
	}
	/**/
	return positions;
}

