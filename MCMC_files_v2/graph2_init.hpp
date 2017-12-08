#ifndef GRAPH2_HPP
#define GRAPH2_HPP

#include <iostream>
#include <string>
#include "mvector.hpp"
#include "graphs.hpp"
#include "graph1.hpp"
//#include "graph3.hpp"

class Graph2
{
	public:
	Graph2(size_t n, size_t q) :
		n(n), q(q), J(xstd::mshape<4>(n, n, q, q)), h(xstd::mshape<2>(n, q)) {}

	Graph2(Graph1 const & g1);

	size_t n, q;
	xstd::mvector<4, double> J;
	xstd::mvector<2, double> h;

	//void randomize(double beta = 1.);

	//void randomize_gauss(double beta = 1.);

	void read(std::istream & is);

	std::ostream & print_distribution(std::ostream & os);

	std::ostream & sample_from_distribution(std::ostream & os, size_t m);

	std::ostream & sample_from_distribution_montecarlo(std::ostream & os, size_t m, size_t mc_iters0, size_t mc_iters, std::string const & out_energies_name, long int seed);

	std::ostream & sample_from_distribution_montecarlo_init(std::ostream & os, size_t m, size_t mc_iters0, size_t mc_iters, std::string const & out_energies_name, int *initial_conf, double *tot_de_record, double *tot_de_record2);

	std::ostream & print_parameters(std::ostream & os);

	void print_parameters(FILE * of);
};

#endif // GRAPH2_HPP
