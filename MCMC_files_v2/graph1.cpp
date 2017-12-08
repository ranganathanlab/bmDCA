#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdio>
#include "mvector.hpp"
//#include "gasdev.hpp"
#include "graphs.hpp"
#include "graph1.hpp"
#include "graph2_init.hpp"
//#include "graph3.hpp"

using namespace std;
using namespace xstd;

extern ostream & log_out;

Graph1::Graph1(Graph2 const & g2) : n(g2.n), q(g2.q), K(mshape<4>(n, n, q, q))
{
	log_out << "converting Ks" << endl;
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					K[i][j][yi][yj] = g2.J[i][j][yi][yj] + 1. / (n - 1) * (g2.h[i][yi] + g2.h[j][yj]);
					K[j][i][yj][yi] = K[i][j][yi][yj];
				}
			}
		}
	}
	log_out << "conversion done" << endl;
}

Graph1::Graph1(Graph3 const & g3) : n(g3.n), q(g3.q), K(mshape<4>(n, n, q, q))
{
	log_out << "converting Ks" << endl;
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					K[i][j][yi][yj] = g3.L[i][j][yi][yj] + 1. / (n - 1) * (g3.l[i][yi] + g3.l[j][yj]);
					K[j][i][yj][yi] = K[i][j][yi][yj];
				}
			}
		}
	}
	log_out << "conversion done" << endl;
}

void Graph1::randomize(double beta)
{
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					K[i][j][yi][yj] = beta * (2 * drand48() - 1);
					K[j][i][yj][yi] = K[i][j][yi][yj];
				}
			}
		}
	}
}

void Graph1::randomize_gauss(double beta)
{
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			log_out << "\rrandomizing " << i << "," << j << "          ";
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					K[i][j][yi][yj] = beta * gasdev();
					K[j][i][yj][yi] = K[i][j][yi][yj];
				}
			}
		}
	}
	log_out << endl;
}

void Graph1::randomize_gauss2(double beta, double het, double base)
{
	assert(het >= 0 and het < 1.);
	vector<double> bi(n);
	vector<size_t> qi(n);
	for (size_t i = 0; i < n; ++i) {
		bi[i] = sqrt(beta) * ((1. - het) + 2 * het * drand48());
		qi[i] = size_t(q * drand48());
	}

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			log_out << "\rrandomizing " << i << "," << j << "          ";
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					K[i][j][yi][yj] = bi[i] * bi[j] * (base * ((qi[i] == yi) + (qi[j] == yj)) + gasdev());
					K[j][i][yj][yi] = K[i][j][yi][yj];
				}
			}
		}
	}
	log_out << endl;
}

ostream & Graph1::print_distribution(ostream & os)
{
	vector<size_t> conf(n);

	double norm = 0;
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += K[i][j][conf[i]][conf[j]];
			}
		}
		norm += exp(x);
		size_t j = 0;
		while (j < n && ++conf[j] == q) {
			conf[j] = 0;
			j++;
		}
		if (j == n) {
			break;
		}
	}
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += K[i][j][conf[i]][conf[j]];
			}
		}
		os << "G1 " << exp(x) / norm << endl;
		size_t j = 0;
		while (j < n && ++conf[j] == q) {
			conf[j] = 0;
			j++;
		}
		if (j == n) {
			break;
		}
	}
	return os;
}

ostream & Graph1::print_parameters(ostream & os)
{
	log_out << "printing parameters K" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					os << "K " << i << " " << j << " " << yi << " " << yj << " " << K[i][j][yi][yj] << endl;
				}
			}
		}
	}
	log_out << "done" << endl;
	return os;
}

void Graph1::print_parameters(FILE * of)
{
	log_out << "printing parameters K" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					fprintf(of, "K %lu %lu %lu %lu %g\n", i, j, yi, yj, K[i][j][yi][yj]);
				}
			}
		}
	}
	log_out << "done" << endl;
}
