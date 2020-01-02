#include <armadillo>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "graph.hpp"
#include "mvector.hpp"

using namespace std;
using namespace xstd;

std::ostream& log_out = std::cerr;

void
Graph::load(potts_model model)
{
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      vector<double> sumJi(q);
      vector<double> sumJj(q);
      for (size_t yi = 0; yi < q; yi++) {
        for (size_t yj = 0; yj < q; yj++) {
          J[i][j][yi][yj] = model.J(i, j)(yi, yj);
          J[j][i][yj][yi] = J[i][j][yi][yj];
          sumJi[yi] += J[i][j][yi][yj];
          sumJi[yj] += J[i][j][yi][yj];
        }
      }
    }
  }
  for (size_t i = 0; i < n; ++i) {
    double sumh = 0;
    for (size_t yi = 0; yi < q; ++yi) {
      h[i][yi] = model.h(yi, i);
      sumh += h[i][yi];
    }
  }
}

void
Graph::read(istream& is)
{
  log_out << "reading Graph" << endl;
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
  }
  log_out << "done" << endl;
}

ostream&
Graph::print_distribution(ostream& os)
{
  vector<size_t> conf(n);

  double norm = 0;
  while (true) {
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

ostream&
Graph::sample_distribution(ostream& os, size_t m)
{
  vector<size_t> conf(n);

  size_t num = size_t(round(pow(double(q), double(n))));

  vector<double> cumulative(num + 1);
  std::cout << "flag 1" << std::endl;

  size_t c = 1;
  while (true) {
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

    size_t j = 0;
    while (j < n && ++conf[j] == q) {
      conf[j] = 0;
      j++;
    }
    if (j == n) {
      break;
    } else {
    }
  }
  assert(c == num + 1);

  double norm = cumulative[num];

  size_t s = 0;
  while (s < m) {
    double x = norm * drand48();

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

void
Graph::sample_mcmc(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   long int seed)
{
  srand48(seed);

  size_t ts = 0;
  vector<size_t> conf(n);
  for (size_t i = 0; i < n; ++i) {
    conf[i] = size_t(q * drand48());
    assert(conf[i] < q);
  }

  double en = 0.;
  for (size_t i = 0; i < n; ++i) {
    en -= h[i][conf[i]];
    for (size_t j = i + 1; j < n; ++j) {
      en -= J[i][j][conf[i]][conf[j]];
    }
  }

  double tot_de = 0;
  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * drand48());
    size_t dq = 1 + size_t((q - 1) * drand48());

    size_t q0 = conf[i];
    size_t q1 = (q0 + dq) % q;

    double e0 = -h[i][q0];
    for (size_t j = 0; j < n; ++j)
      if (j != i) {
        e0 -= J[i][j][q0][conf[j]];
      }
    double e1 = -h[i][q1];
    for (size_t j = 0; j < n; ++j)
      if (j != i) {
        e1 -= J[i][j][q1][conf[j]];
      }
    double de = e1 - e0;
    if ((de < 0) || (drand48() < exp(-de))) {
      conf[i] = q1;
      tot_de += de;
    }
  }
  en += tot_de;
  tot_de = 0.;
  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * drand48());
      size_t dq = 1 + size_t((q - 1) * drand48());

      size_t q0 = conf[i];
      size_t q1 = (q0 + dq) % q;

      double e0 = -h[i][q0];
      for (size_t j = 0; j < n; ++j)
        if (j != i) {
          e0 -= J[i][j][q0][conf[j]];
        }
      double e1 = -h[i][q1];
      for (size_t j = 0; j < n; ++j)
        if (j != i) {
          e1 -= J[i][j][q1][conf[j]];
        }
      double de = e1 - e0;
      if ((de < 0) || (drand48() < exp(-de))) {
        conf[i] = q1;
        tot_de += de;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      (*ptr)(s, i) = conf[i];
    }
  }
  log_out << "sampled " << m << " [de=" << tot_de << "]" << endl;
  return;
}

ostream&
Graph::sample_mcmc(ostream& os,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   string const& out_energies_name,
                   long int seed)
{
  srand48(seed);

  size_t ts = 0;
  vector<size_t> conf(n);
  for (size_t i = 0; i < n; ++i) {
    conf[i] = size_t(q * drand48());
    assert(conf[i] < q);
  }

  log_out << "computing initial energy... ";
  double en = 0.;
  for (size_t i = 0; i < n; ++i) {
    en -= h[i][conf[i]];
    for (size_t j = i + 1; j < n; ++j) {
      en -= J[i][j][conf[i]][conf[j]];
    }
  }
  log_out << "done." << endl;

  log_out << "initialize montecarlo sampling... ";
  double tot_de = 0;
  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * drand48());
    size_t dq = 1 + size_t((q - 1) * drand48());

    size_t q0 = conf[i];
    size_t q1 = (q0 + dq) % q;

    double e0 = -h[i][q0];
    for (size_t j = 0; j < n; ++j)
      if (j != i) {
        e0 -= J[i][j][q0][conf[j]];
      }
    double e1 = -h[i][q1];
    for (size_t j = 0; j < n; ++j)
      if (j != i) {
        e1 -= J[i][j][q1][conf[j]];
      }
    double de = e1 - e0;
    if ((de < 0) || (drand48() < exp(-de))) {
      conf[i] = q1;
      tot_de += de;
    }
  }
  log_out << " [tot_de=" << tot_de << "] done." << endl;
  en += tot_de;
  tot_de = 0.;
  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * drand48());
      size_t dq = 1 + size_t((q - 1) * drand48());

      size_t q0 = conf[i];
      size_t q1 = (q0 + dq) % q;

      double e0 = -h[i][q0];
      for (size_t j = 0; j < n; ++j)
        if (j != i) {
          e0 -= J[i][j][q0][conf[j]];
        }
      double e1 = -h[i][q1];
      for (size_t j = 0; j < n; ++j)
        if (j != i) {
          e1 -= J[i][j][q1][conf[j]];
        }
      double de = e1 - e0;
      if ((de < 0) || (drand48() < exp(-de))) {
        conf[i] = q1;
        tot_de += de;
      }
    }
    log_out << "\rs=" << ++ts << "/" << m << " de=" << tot_de
            << "                 ";
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
  return os;
}

ostream&
Graph::initialize_mcmc(ostream& os,
                       size_t m,
                       size_t mc_iters0,
                       size_t mc_iters,
                       string const& out_energies_name,
                       int* initial_conf,
                       double* tot_de_record,
                       double* tot_de_record2)
{

  srand48(time(NULL));

  size_t ts = 0;
  vector<size_t> conf(n);
  for (size_t i = 0; i < n; ++i) {
    conf[i] = initial_conf[i];
    assert(conf[i] < q);
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
    for (size_t j = 0; j < n; ++j)
      if (j != i) {
        e0 -= J[i][j][q0][conf[j]];
      }
    double e1 = -h[i][q1];
    for (size_t j = 0; j < n; ++j)
      if (j != i) {
        e1 -= J[i][j][q1][conf[j]];
      }
    double de = e1 - e0;
    if ((de < 0) || (drand48() < exp(-de))) {
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
      for (size_t j = 0; j < n; ++j)
        if (j != i) {
          e0 -= J[i][j][q0][conf[j]];
        }
      double e1 = -h[i][q1];
      for (size_t j = 0; j < n; ++j)
        if (j != i) {
          e1 -= J[i][j][q1][conf[j]];
        }
      double de = e1 - e0;
      if ((de < 0) || (drand48() < exp(-de))) {
        // cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
        conf[i] = q1;
        tot_de += de;
      }
      tot_de_record[k] += tot_de;
      tot_de_record2[k] += tot_de * tot_de;
    }
    log_out << "\rs=" << ++ts << "/" << m << " de=" << tot_de
            << "                 ";
    out_energies << s + 1 << " " << en + tot_de << endl;
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

ostream&
Graph::print_parameters(ostream& os)
{
  log_out << "printing parameters J" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      for (size_t yi = 0; yi < q; ++yi) {
        for (size_t yj = 0; yj < q; ++yj) {
          os << "J " << i << " " << j << " " << yi << " " << yj << " "
             << J[i][j][yi][yj] << endl;
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

void
Graph::print_parameters(FILE* of)
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
