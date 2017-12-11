#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <ctime>
#include "mvector.hpp"
//#include "gasdev.hpp"
#include "graph1.hpp"
//#include "graph2.hpp"
#include "graph2_init.hpp"
//#include "graph3.hpp"

using namespace std;
using namespace xstd;

extern ostream & log_out;

Graph2::Graph2(Graph1 const & g1) : n(g1.n), q(g1.q), J(mshape<4>(n, n, q, q)), h(mshape<2>(n, q))
{
	log_out << "converting Js" << endl;
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			/*cerr << "K:" << endl;
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					cerr << g1.K[i][j][yi][yj] << " ";
				}
				cerr << endl;
			}*/
			vector<double> sumKi(q), sumKj(q);
			double sumKall = 0.;
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					sumKi[yi] += g1.K[i][j][yi][yj];
					sumKj[yj] += g1.K[i][j][yi][yj];
					sumKall += g1.K[i][j][yi][yj];
				}
			}
			/*cerr << endl;
			cerr << "sumKi: "; for (size_t yi = 0; yi < q; ++yi) { cerr << sumKi[yi] << " "; } cerr << endl;
			cerr << "sumKj: "; for (size_t yj = 0; yj < q; ++yj) { cerr << sumKj[yj] << " "; } cerr << endl;
			cerr << "sumKall: " << sumKall << endl;*/
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					J[i][j][yi][yj] = g1.K[i][j][yi][yj] - (1. / q) * (sumKi[yi] + sumKj[yj]) + 1. / (q * q) * sumKall; 
					//J[i][j][yi][yj] /= 5; // XXX
					J[j][i][yj][yi] = J[i][j][yi][yj];
				}
			}
			/*cerr << endl;
			cerr << "J:" << endl;
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					cerr << J[i][j][yi][yj] << " ";
				}
				cerr << endl;
			}
			cerr << endl;*/
			vector<double> sumJi(q), sumJj(q);
			double sumJall = 0.;
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					sumJi[yi] += J[i][j][yi][yj];
					sumJj[yj] += J[i][j][yi][yj];
					sumJall += J[i][j][yi][yj];
				}
			}
			/*cerr << endl;
			cerr << "sumJi: "; for (size_t yi = 0; yi < q; ++yi) { cerr << sumJi[yi] << " "; } cerr << endl;
			cerr << "sumJj: "; for (size_t yj = 0; yj < q; ++yj) { cerr << sumJj[yj] << " "; } cerr << endl;
			cerr << "sumJall: " << sumJall << endl;*/
			for (size_t y = 0; y < q; ++y) {
				assert(fabs(sumJi[y]) < 1e-10);
				assert(fabs(sumJj[y]) < 1e-10);
			}
			//cout<<sumJall<<endl;
			assert(fabs(sumJall) < 1e-10);
		}
	}
	log_out << "converting Hs" << endl;
#pragma omp parallel for
	for (size_t i = 0; i < n; ++i) {
		vector<double> sumKi(q);
		double sumKall = 0.;
		for (size_t yi = 0; yi < q; ++yi) {
			for (size_t j = 0; j < n; ++j) if (j != i) {
				for (size_t yj = 0; yj < q; ++yj) {
					sumKi[yi] += g1.K[i][j][yi][yj];
					sumKall += g1.K[i][j][yi][yj];
				}
			}
		}
		for (size_t yi = 0; yi < q; ++yi) {
			h[i][yi] = 1. / q * sumKi[yi] - 1. / (q * q) * sumKall;
		}
		double sumhall = 0.;
		for (size_t yi = 0; yi < q; ++yi) {
			sumhall += h[i][yi];
		}
		//cerr << "sumhall=" << sumhall << endl;
		assert(fabs(sumhall) < 1e-10);
	}
	log_out << "conversion done" << endl;
}

void Graph2::read(istream & is)
{
	log_out << "reading Graph2" << endl;
	string tok;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			vector<double> sumJi(q);
			vector<double> sumJj(q);
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					is >> tok;
					assert(tok == "J");
					is >> tok;
					assert(atoi(tok.c_str()) == int(i));
					is >> tok;
					assert(atoi(tok.c_str()) == int(j));
					is >> tok;
					assert(atoi(tok.c_str()) == int(yi));
					is >> tok;
					assert(atoi(tok.c_str()) == int(yj));
					is >> tok;
					J[i][j][yi][yj] = atof(tok.c_str());
					J[j][i][yj][yi] = J[i][j][yi][yj];
					sumJi[yi] += J[i][j][yi][yj];
					sumJi[yj] += J[i][j][yi][yj];
				}
			}
			for (size_t y = 0; y < q; ++y) {
				//cerr << i << " " << j << " " << " sumJi[y]=" << sumJi[y] << endl;
				//cerr << i << " " << j << " " << " sumJj[y]=" << sumJj[y] << endl;
				assert(fabs(sumJi[y]) < 1e10);
				assert(fabs(sumJj[y]) < 1e10);
			}
		}
	}
	for (size_t i = 0; i < n; ++i) {
		double sumh = 0;
		for (size_t yi = 0; yi < q; ++yi) {
			is >> tok;
			assert(tok == "h");
			is >> tok;
			assert(atoi(tok.c_str()) == int(i));
			is >> tok;
			assert(atoi(tok.c_str()) == int(yi));
			is >> tok;
			h[i][yi] = atof(tok.c_str());
			sumh += h[i][yi];
		}
		//assert(fabs(sumh) < 1e-3);
	}
	log_out << "done" << endl;
}

ostream & Graph2::print_distribution(ostream & os)
{
	vector<size_t> conf(n);

	double norm = 0;
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			x += h[i][conf[i]];
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += J[i][j][conf[i]][conf[j]];
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
			x += h[i][conf[i]];
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += J[i][j][conf[i]][conf[j]];
			}
		}
		os << "G2 " << exp(x) / norm << endl;
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


ostream & Graph2::sample_from_distribution(ostream & os, size_t m)
{
	vector<size_t> conf(n);

	size_t num = size_t(round(pow(double(q), double(n))));

	vector<double> cumulative(num + 1);

	size_t c = 1;
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			x += h[i][conf[i]];
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += J[i][j][conf[i]][conf[j]];
			}
		}
		double nnp = exp(x);
		cumulative[c] = cumulative[c - 1];
		cumulative[c++] += nnp;
		//cout << c - 2 << ": ( "; for (size_t i = 0; i < n; ++i) cout << conf[i] << " "; cout << ") : " << exp(x) << " " << cumulative[c - 1] << endl;

		size_t j = 0;
		while (j < n && ++conf[j] == q) {
			conf[j] = 0;
			j++;
		}
		if (j == n) {
			break;
		}
	}
	assert(c == num + 1);

	double norm = cumulative[num];

	size_t s = 0;
	while (s < m) {
		double x = norm * drand48();
		//double x = norm * (double(s) / m);

		size_t i0 = 0;
		size_t i1 = num;
		while (i1 != i0 + 1) {
			size_t i = i0 + (i1 - i0) / 2;
			if (x < cumulative[i]) {
				i1 = i;
			} else {
				i0 = i;
			}
		}
		//cout << "x=" << x << " i0=" << i0 << " : ";
		size_t qq = i0;
		for (size_t i = 0; i < n; ++i) {
			os << qq % q;
			if (i < n - 1) {
				os << " ";
			} else {
				os << endl;
			}
			qq /= q;
		}
		++s;
	}
	return os;
}

ostream & Graph2::sample_from_distribution_montecarlo(ostream & os, size_t m, size_t mc_iters0, size_t mc_iters, string const & out_energies_name, long int seed)
{

        srand48(seed);
//        srand48(time(NULL));

	size_t ts = 0;
	vector<size_t> conf(n);
	for (size_t i = 0; i < n; ++i) {
		conf[i] = size_t(q * drand48());
		assert(conf[i] < q);
		//cerr << conf[i] << " ";
	}
	ofstream out_energies(out_energies_name.c_str());
	log_out << "computing initial energy... ";
	double en = 0.;
	for (size_t i = 0; i < n; ++i) {
		en -= h[i][conf[i]];
		for (size_t j = i + 1; j < n; ++j) {
			en -= J[i][j][conf[i]][conf[j]];
		}
	}
	log_out << "done." << endl;
	out_energies << "-1 " << en << endl;

	log_out << "initialize montecarlo sampling... ";
	double tot_de = 0;
	for (size_t k = 0; k < mc_iters0; ++k) {
		size_t i = size_t(n * drand48());
		size_t dq = 1 + size_t((q - 1) * drand48());

		size_t q0 = conf[i];
		size_t q1 = (q0 + dq) % q;

		double e0 = -h[i][q0];
		for (size_t j = 0; j < n; ++j) if (j != i) {
			e0 -= J[i][j][q0][conf[j]];
		}
		double e1 = -h[i][q1];
		for (size_t j = 0; j < n; ++j) if (j != i) {
			e1 -= J[i][j][q1][conf[j]];
		}
		double de = e1 - e0;
		if ((de < 0) || (drand48() < exp(-de))) {
			//cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
			conf[i] = q1;
			tot_de += de;
		}
	}
	log_out << " [tot_de=" << tot_de << "] done." << endl;
	en += tot_de;
	out_energies << "0 " << en << endl;
	tot_de = 0.;
	for (size_t s = 0; s < m; ++s) {
		for (size_t k = 0; k < mc_iters; ++k) {
			size_t i = size_t(n * drand48());
			size_t dq = 1 + size_t((q - 1) * drand48());

			size_t q0 = conf[i];
			size_t q1 = (q0 + dq) % q;

			double e0 = -h[i][q0];
			for (size_t j = 0; j < n; ++j) if (j != i) {
				e0 -= J[i][j][q0][conf[j]];
			}
			double e1 = -h[i][q1];
			for (size_t j = 0; j < n; ++j) if (j != i) {
				e1 -= J[i][j][q1][conf[j]];
			}
			double de = e1 - e0;
			if ((de < 0) || (drand48() < exp(-de))) {
				//cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
				conf[i] = q1;
				tot_de += de;
			}
		}
		log_out << "\rs=" << ++ts << "/" << m << " de=" << tot_de << "                 ";
		out_energies << s+1 << " " << en + tot_de << endl;
		for (size_t i = 0; i < n; ++i) {
			os << conf[i];
			if (i < n - 1) {
				os << " ";
			} else {
				os << endl;
			}
		}
	}
	log_out << endl;
	out_energies.close();
	return os;
}



ostream & Graph2::sample_from_distribution_montecarlo_init(ostream & os, size_t m, size_t mc_iters0, size_t mc_iters, string const & out_energies_name, int * initial_conf, double *tot_de_record,double *tot_de_record2)
{

        srand48(time(NULL));

	size_t ts = 0;
	vector<size_t> conf(n);
	for (size_t i = 0; i < n; ++i) {
		conf[i] = initial_conf[i];
		assert(conf[i] < q);
		//cerr << conf[i] << " ";
	}
	ofstream out_energies(out_energies_name.c_str());
	log_out << "computing initial energy... ";
	double en = 0.;
	for (size_t i = 0; i < n; ++i) {
		en -= h[i][conf[i]];
		for (size_t j = i + 1; j < n; ++j) {
			en -= J[i][j][conf[i]][conf[j]];
		}
	}
	log_out << "done." << endl;
	out_energies << "-1 " << en << endl;

	log_out << "initialize montecarlo sampling... ";
	double tot_de = 0;
	for (size_t k = 0; k < mc_iters0; ++k) {
		size_t i = size_t(n * drand48());
		size_t dq = 1 + size_t((q - 1) * drand48());

		size_t q0 = conf[i];
		size_t q1 = (q0 + dq) % q;

		double e0 = -h[i][q0];
		for (size_t j = 0; j < n; ++j) if (j != i) {
			e0 -= J[i][j][q0][conf[j]];
		}
		double e1 = -h[i][q1];
		for (size_t j = 0; j < n; ++j) if (j != i) {
			e1 -= J[i][j][q1][conf[j]];
		}
		double de = e1 - e0;
		if ((de < 0) || (drand48() < exp(-de))) {
			//cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
			conf[i] = q1;
			tot_de += de;
		}
	}
	log_out << " [tot_de=" << tot_de << "] done." << endl;
	en += tot_de;
	out_energies << "0 " << en << endl;
	tot_de = 0.;
	for (size_t s = 0; s < m; ++s) {
printf("b");
fflush(0);
		for (size_t k = 0; k < mc_iters; ++k) {
printf("a");
fflush(0);
			size_t i = size_t(n * drand48());
			size_t dq = 1 + size_t((q - 1) * drand48());

			size_t q0 = conf[i];
			size_t q1 = (q0 + dq) % q;

			double e0 = -h[i][q0];
			for (size_t j = 0; j < n; ++j) if (j != i) {
				e0 -= J[i][j][q0][conf[j]];
			}
			double e1 = -h[i][q1];
			for (size_t j = 0; j < n; ++j) if (j != i) {
				e1 -= J[i][j][q1][conf[j]];
			}
			double de = e1 - e0;
			if ((de < 0) || (drand48() < exp(-de))) {
				//cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
				conf[i] = q1;
				tot_de += de;				
			}
			tot_de_record[k]+=tot_de;
			tot_de_record2[k]+=tot_de*tot_de;
		}
		log_out << "\rs=" << ++ts << "/" << m << " de=" << tot_de << "                 ";
		out_energies << s+1 << " " << en + tot_de << endl;
		for (size_t i = 0; i < n; ++i) {
			os << conf[i];
			if (i < n - 1) {
				os << " ";
			} else {
				os << endl;
			}
		}
	}
printf("c");
fflush(0);
	log_out << endl;
	out_energies.close();
	return os;
}


ostream & Graph2::print_parameters(ostream & os)
{
	log_out << "printing parameters J" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					os << "J " << i << " " << j << " " << yi << " " << yj << " " << J[i][j][yi][yj] << endl;
				}
			}
		}
	}
	log_out << "printing parameters H" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t yi = 0; yi < q; ++yi) {
			os << "h " << i << " " << yi << " " << h[i][yi] << endl;
		}
	}
	log_out << "done" << endl;
	return os;
}

void Graph2::print_parameters(FILE * of)
{
	log_out << "printing parameters J" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					fprintf(of, "J %lu %lu %lu %lu %g\n", i, j, yi, yj, J[i][j][yi][yj]);
				}
			}
		}
	}
	log_out << "printing parameters H" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t yi = 0; yi < q; ++yi) {
			fprintf(of, "h %lu %lu %g\n", i, yi, h[i][yi]);
		}
	}
	log_out << "done" << endl;
}
