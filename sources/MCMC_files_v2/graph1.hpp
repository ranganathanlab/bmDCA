#ifndef GRAPH1_HPP
#define GRAPH1_HPP

#include <cassert>
#include <iostream>
#include "mvector.hpp"
#include "graphs.hpp"
#include "graph2_init.hpp"
//#include "graph3.hpp"

class Graph1
{
	public:
	Graph1(size_t n, size_t q) :
		n(n), q(q), K(xstd::mshape<4>(n, n, q, q)) {}

	Graph1(Graph2 const & g2);

	Graph1(Graph3 const & g3);

	size_t n, q;
	xstd::mvector<4, double> K;

	void randomize(double beta = 1.);

	void randomize_gauss(double beta = 1.);

	void randomize_gauss2(double beta = 1., double het = 0.5, double base = 0.);

	std::ostream & print_distribution(std::ostream & os);

	std::ostream & print_parameters(std::ostream & os);

	void print_parameters(FILE * of);
};

#endif // GRAPH1_HPP
