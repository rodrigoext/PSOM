/*
 * Watershed.cpp
 *
 *  Created on: 15/05/2016
 *      Author: rodrigo
 */

#include <Model/Utils/Watershed.h>
#include <Model/Utils/IO.h>
#include <queue>
#include <iostream>
// Níveis de discretização da imagem
#define LEVELS 31
#define MASK -2 /* initial value of a threshold level */
#define WSHED 0 /* value of the pixels belonging to the watersheds */
#define INIT -1 /* initial value of im, */

Watershed::Watershed() {
	// TODO Auto-generated constructor stub

}

Watershed::~Watershed() {
	// TODO Auto-generated destructor stub
}

Eigen::MatrixXf Watershed::transform(const Eigen::MatrixXf &input) {
    Eigen::MatrixXf imi = input*20;

	std::cout << "input" << std::endl;
	std::cout << imi << std::endl;

	int rows = imi.rows();
	int cols = imi.cols();
	int units = rows * cols;
	IO * io;
	io = new IO();
	//imi = imi * LEVELS;
	//io->MapMinMax(imi);
	for (int l = 0; l < rows; l++) {
		for (int c = 0; c < cols; c++) {
			// imi(l, c) = ((int) (imi(l, c) + 0.5)) + 1;
			imi(l, c) = (int) ((imi(l, c) + 1));
		}
	}
	std::cout << imi.rows() << " || " << imi.cols() << std::endl;
	std::cout << imi << std::endl;
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
        std::cout << " (" << row_min << "|" << col_min << ") ";
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
    std::cout << std::endl;

    //std::cout << "temp_mat" << std::endl;
    //std::cout << temp_mat << std::endl;

    std::cout << "sorted: " << hmin << " | " << hmax  << std::endl;

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
            if (last_pixel > (units-1))
				break;
		}
        //std::cout << "imd" << std::endl;
        //std::cout << imd << std::endl;
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
            std::cout << "second while" << std::endl;
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
	std::cout << "labels" << std::endl;
	std::cout << imo << std::endl;
	return imo;
}

Eigen::MatrixXf Watershed::transform_v2(const Eigen::MatrixXf & input) {
    //TODO: verifica se precisa transformar para escala cinza
    Eigen::MatrixXf data(input*18);
    Eigen::MatrixXf lab(data.rows(), data.cols());
    Eigen::MatrixXf dist(data.rows(), data.cols());
    int units = data.cols() * data.cols();
    int curlab = 0;

    int hmin = (int)data.minCoeff();
    int hmax = (int)data.maxCoeff();

    for (int i = 0; i < data.rows(); ++i)
        for (int j = 0; j < data.cols(); ++j) {
            lab(i, j) = INIT;
            dist(i, j) = 0;
        }

    // SORT pixels in increasing order of grey values
    std::vector<pos> position_sorted;
    Eigen::VectorXi pixels_sorted(data.rows() * data.cols());
    Eigen::MatrixXf temp_input(data);
    std::cout << temp_input << std::endl;
    float temp_min;
    for (int i = 0; i < units; ++i) {
        temp_min = temp_input(0, 0);
        pos p;
        p.row = 0;
        p.col = 0;
        for (int r = 0; r < data.rows(); ++r) {
            for (int c = 0; c < data.cols(); ++c) {
                if (temp_min > temp_input(r, c)) {
                    temp_min = temp_input(r, c);
                    p.row = r;
                    p.col = c;
                }
            }
        }
        position_sorted.push_back(p);
        pixels_sorted[i] = (int)temp_min;
        temp_input(p.row, p.col) = hmax + 1;
    }

    // Start Flooding
    std::queue<pos> fifo;
    int last_pixel = 0;
    pos p;
    for (int h = hmin; h < hmax; ++h) {
        // mask all pixels at level h
        while (h == pixels_sorted(last_pixel)) {
            p = position_sorted.at(last_pixel);
            lab(p.row, p.col) = MASK;
            std::vector<pos> pixel_neighbour = vizinhanca4(p.row, p.col, data.rows(), data.cols());
            for (auto pn : pixel_neighbour) {
                if(lab(pn.row, pn.col) > 0 || lab(pn.row, pn.col) == WSHED) {
                    dist(p.row, p.col) = 1;
                    fifo.push(pn);
                }
            }
            last_pixel++;
            if (last_pixel > (pixels_sorted.size() - 1))
                break;
        }

        // extend basians
        int curdist = 1;
        pos fictitious;
        fictitious.row = -1;
        fictitious.col = -1;
        fifo.push(fictitious);

        while (1) {
            p = fifo.front();
            fifo.pop();
            if (p.row == fictitious.row && p.col == fictitious.col) {
                if (fifo.empty()) {
                    break;
                } else {
                    fifo.push(fictitious);
                    curdist++;
                    p = fifo.front();
                    fifo.pop();
                }
            }
            //labelling p by inspecting neighbours
            std::vector<pos> pixel_neighbour = vizinhanca4(p.row, p.col, data.rows(), data.cols());
            for (auto q : pixel_neighbour) {
                if (dist(q.row, q.col) < curdist && (lab(q.row, q.col) > 0 || lab(q.row, q.col) == WSHED)) {
                    //q belongs to an existing basin or to watersheds
                    if (lab(q.row, q.col) > 0) {
                        if (lab(p.row, p.col) == MASK || lab(p.row, p.col) == WSHED) {
                            std::cout << "lab atrib" << std::endl;
                            lab(p.row, p.col) = lab(q.row, q.col);
                        } else if (lab(p.row, p.col) != lab(q.row, q.col)) {
                            std::cout << "lab atrib" << std::endl;
                            lab(p.row, p.col) = WSHED;
                        }
                    } else if (lab(p.row, p.col) == MASK) {
                        std::cout << "lab atrib" << std::endl;
                        lab(p.row, p.col) = WSHED;
                    }
                } else if (lab(q.row, q.col) == MASK && dist(q.row, q.col) == 0) {
                    dist(q.row, q.col) = curdist + 1;
                    fifo.push(q);
                }
            }
        }

        // detect and process new minima at level h
        int last_pixel_temp = 0;
        pos p2;
        while (h == pixels_sorted(last_pixel_temp)) {
            p2 = position_sorted.at(last_pixel_temp);
            dist(p2.row, p2.col) = 0;
            if (lab(p2.row, p2.col) == MASK) {
                curlab = curlab + 1;
                fifo.push(p2);
                std::cout << "lab atrib" << std::endl;
                lab(p2.row, p2.col) = curlab;
                while (!fifo.empty()) {
                    pos q = fifo.front();
                    fifo.pop();
                    std::vector<pos> pixel_neighbour = vizinhanca4(q.row, q.col, data.rows(), data.cols());
                    for (auto r : pixel_neighbour) {
                        if (lab(r.row, r.col) == MASK) {
                            fifo.push(r);
                            lab(r.row, r.col) = curlab;
                        }
                    }
                }
            }
            last_pixel_temp++;
            if (last_pixel_temp > (pixels_sorted.size() - 1))
                break;
        }
    }

    std::cout << "v2 end" << std::endl;
    std::cout << lab << std::endl;
    return lab;
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
