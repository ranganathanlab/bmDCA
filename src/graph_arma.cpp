#include <armadillo>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>

#include "graph_arma.hpp"
#include "pcg_random.hpp"

using namespace std;

std::ostream& log_out = std::cout;

void
Graph::load(potts_model* model) {
  params = model;
}

ostream&
Graph::print_distribution(ostream& os)
{
  vector<size_t> conf(n);

  double norm = 0;
  while (true) {
    double x = 0;
    for (size_t i = 0; i < n; ++i) {
      x += params->h.at(conf[i], i);
    }
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        x += params->J.at(i, j).at(conf[i], conf[j]);
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
      x += params->h.at(conf[i], i);
    }
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        x += params->J.at(i, j).at(conf[i], conf[j]);
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
};

// ostream&
// Graph::sample_distribution(ostream& os, size_t m)
// {
//   vector<size_t> conf(n);
//
//   size_t num = size_t(round(pow(double(q), double(n))));
//
//   vector<double> cumulative(num + 1);
//   std::cout << "flag 1" << std::endl;
//
//   size_t c = 1;
//   while (true) {
//     double x = 0;
//     for (size_t i = 0; i < n; ++i) {
//       x += h[i][conf[i]];
//     }
//     for (size_t i = 0; i < n; ++i) {
//       for (size_t j = i + 1; j < n; ++j) {
//         x += J[i][j][conf[i]][conf[j]];
//       }
//     }
//     double nnp = exp(x);
//     cumulative[c] = cumulative[c - 1];
//     cumulative[c++] += nnp;
//
//     size_t j = 0;
//     while (j < n && ++conf[j] == q) {
//       conf[j] = 0;
//       j++;
//     }
//     if (j == n) {
//       break;
//     } else {
//     }
//   }
//   assert(c == num + 1);
//
//   double norm = cumulative[num];
//
//   size_t s = 0;
//   while (s < m) {
//     double x = norm * drand48();
//
//     size_t i0 = 0;
//     size_t i1 = num;
//     while (i1 != i0 + 1) {
//       size_t i = i0 + (i1 - i0) / 2;
//       if (x < cumulative[i]) {
//         i1 = i;
//       } else {
//         i0 = i;
//       }
//     }
//     size_t qq = i0;
//     for (size_t i = 0; i < n; ++i) {
//       os << qq % q;
//       if (i < n - 1) {
//         os << " ";
//       } else {
//         os << endl;
//       }
//       qq /= q;
//     }
//     ++s;
//   }
//   return os;
// }

void
Graph::sample_mcmc(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   long int seed,
                   double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  // size_t ts = 0;
  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf.at(i) = size_t(q * uniform(rng));
    assert(conf.at(i) < q);
  }

  double en = 0.;
  for (size_t i = 0; i < n; ++i) {
    en -= params->h.at(conf.at(i), i);
    for (size_t j = i + 1; j < n; ++j) {
      en -= params->J.at(i, j).at(conf.at(i), conf.at(j));
    }
  }

  double tot_de = 0;
  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * uniform(rng));
    size_t dq = 1 + size_t((q - 1) * uniform(rng));

    size_t q0 = conf.at(i);
    size_t q1 = (q0 + dq) % q;

    double e0 = -params->h.at(q0, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J.at(j, i).at(conf.at(j), q0);
      } else if (i < j) {
        e0 -= params->J.at(i, j).at(q0, conf.at(j));
      }
    }
    double e1 = -params->h.at(q1, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e1 -= params->J.at(j, i).at(conf.at(j), q1);
      } else if (i < j) {
        e1 -= params->J.at(i, j).at(q1, conf.at(j));
      }
    }
    double de = e1 - e0;
    if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
      conf.at(i) = q1;
      tot_de += de;
    }
  }

  en += tot_de;
  tot_de = 0.;
  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * uniform(rng));
      size_t dq = 1 + size_t((q - 1) * uniform(rng));

      size_t q0 = conf.at(i);
      size_t q1 = (q0 + dq) % q;

      double e0 = -params->h.at(q0, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e0 -= params->J.at(j, i).at(conf.at(j), q0);
        } else if (i < j) {
          e0 -= params->J.at(i, j).at(q0, conf.at(j));
        }
      }
      double e1 = -params->h.at(q1, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e1 -= params->J.at(j, i).at(conf.at(j), q1);
        } else if (i < j) {
          e1 -= params->J.at(i, j).at(q1, conf.at(j));
        }
      }
      double de = e1 - e0;
      if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
        conf.at(i) = q1;
        tot_de += de;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      (*ptr).at(s, i) = conf.at(i);
    }
  }
  // std::string output_string =
  //   "sampled " + std::to_string(m) + " [de=" + std::to_string(en + tot_de) + "]\n";
  // log_out << output_string;
  return;
};

void
Graph::sample_mcmc_init(arma::Mat<int>* ptr,
                        size_t m,
                        size_t mc_iters0,
                        size_t mc_iters,
                        arma::Col<int>* init_ptr,
                        long int seed,
                        double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  // size_t ts = 0;
  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf.at(i) = (*ptr).at(i);
    assert(conf.at(i) < q);
  }

  double en = 0.;
  for (size_t i = 0; i < n; ++i) {
    en -= params->h.at(conf.at(i), i);
    for (size_t j = i + 1; j < n; ++j) {
      en -= params->J.at(i, j).at(conf.at(i), conf.at(j));
    }
  }

  double tot_de = 0;
  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * uniform(rng));
    size_t dq = 1 + size_t((q - 1) * uniform(rng));

    size_t q0 = conf.at(i);
    size_t q1 = (q0 + dq) % q;

    double e0 = -params->h.at(q0, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J.at(j, i).at(conf.at(j), q0);
      } else if (i < j) {
        e0 -= params->J.at(i, j).at(q0, conf.at(j));
      }
    }
    double e1 = -params->h.at(q1, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e1 -= params->J.at(j, i).at(conf.at(j), q1);
      } else if (i < j) {
        e1 -= params->J.at(i, j).at(q1, conf.at(j));
      }
    }
    double de = e1 - e0;
    if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
      conf.at(i) = q1;
      tot_de += de;
    }
  }
  en += tot_de;
  tot_de = 0.;
  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * uniform(rng));
      size_t dq = 1 + size_t((q - 1) * uniform(rng));

      size_t q0 = conf.at(i);
      size_t q1 = (q0 + dq) % q;

      double e0 = -params->h.at(q0, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e0 -= params->J.at(j, i).at(conf.at(j), q0);
        } else if (i < j) {
          e0 -= params->J.at(i, j).at(q0, conf.at(j));
        }
      }
      double e1 = -params->h.at(q1, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e1 -= params->J.at(j, i).at(conf.at(j), q1);
        } else if (i < j) {
          e1 -= params->J.at(i, j).at(q1, conf.at(j));
        }
      }
      double de = e1 - e0;
      if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
        conf.at(i) = q1;
        tot_de += de;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      (*ptr).at(s, i) = conf.at(i);
    }
  }
  // std::string output_string =
  //   "sampled " + std::to_string(m) + " [de=" + std::to_string(tot_de) + "]\n";
  // log_out << output_string;
  return;
};

void
Graph::sample_mcmc_zanella(arma::Mat<int>* ptr,
                           size_t m,
                           size_t mc_iters0,
                           size_t mc_iters,
                           long int seed,
                           double temperature,
                           std::string mode)
{

  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  // Generate random initial sequence.
  // size_t ts = 0;
  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf.at(i) = size_t(q * uniform(rng));
    assert(conf.at(i) < q);
  }

  // // Compute energy of random initial sequence.
  // double en_x = 0.;
  // for (size_t i = 0; i < n; ++i) {
  //   en_x -= h[i][conf[i]];
  //   for (size_t j = i + 1; j < n; ++j) {
  //     en_x -= J[i][j][conf[i]][conf[j]];
  //   }
  // }

  arma::Mat<double> de = arma::Mat<double>(n, q, arma::fill::zeros);
  arma::Mat<double> g = arma::Mat<double>(n, q, arma::fill::zeros);
  double lambda = 0.0;

  for (size_t k = 0; k < mc_iters0; ++k) {

    // Compute initial neighborhood.
    if (k < 1) {
      for (size_t i = 0; i < n; i++) {
        size_t q0 = conf.at(i);
        double e0 = -params->h.at(q0, i);
        for (size_t j = 0; j < n; ++j) {
          if (i > j) {
            e0 -= params->J.at(j, i).at(conf.at(j), q0);
          } else if (i < j) {
            e0 -= params->J.at(i, j).at(q0, conf.at(j));
          }
        }

        for (size_t q1 = 0; q1 < q; q1++) {
          if (q0 == q1) {
            de.at(i, q1) = 0.0;
          } else {
            double e1 = -params->h.at(q1, i);
            for (size_t j = 0; j < n; ++j) {
              if (i > j) {
                e1 -= params->J.at(j, i).at(conf.at(j), q1);
              } else if (i < j) {
                e1 -= params->J.at(i, j).at(q1, conf.at(j));
              }
            }
            de.at(i, q1) = e1 - e0;
          }
        }
      }
    }

    if (mode == "sqrt") {
      g = arma::exp(de * -1.0 / 2.0 / temperature);
      lambda = arma::accu(g) - n; // n*exp(0) needs to be subtracted.
    } else if (mode == "tanh") {
      g = .5 + .5 * arma::tanh(de * -1.0 / 2.0 / temperature);
      lambda = arma::accu(g) - .5 * n;
    }
    g = g / lambda;

    double rand = uniform(rng);
    double r_sum = 0.0;
    size_t i = 0;
    size_t q0 = 0;
    size_t q1 = 0;
    bool done = false;
    for (i = 0; i < n; i++) {
      for (q1 = 0; q1 < q; q1++) {
        if (conf.at(i) == q1) {
          continue;
        } else {
          r_sum += g.at(i, q1);
        }
        if (r_sum > rand) {
          done = true;
          break;
        }
      }
      if (done) {
        break;
      }
    }

    double tmp = de.at(i, q1);
    q0 = conf.at(i);
    conf.at(i) = q1;
    for (size_t pos = 0; pos < n; pos++) {
      for (size_t aa = 0; aa < q; aa++) {
        if (pos != i) {
          if (pos < i) {
            de.at(pos, aa) += params->J.at(pos, i).at(conf.at(pos), q1) -
                              params->J.at(pos, i).at(conf.at(pos), q0) -
                              params->J.at(pos, i).at(aa, q1) +
                              params->J.at(pos, i).at(aa, q0);
          } else {
            de.at(pos, aa) += params->J.at(i, pos).at(q1, conf.at(pos)) -
                              params->J.at(i, pos).at(q0, conf.at(pos)) -
                              params->J.at(i, pos).at(q1, aa) +
                              params->J.at(i, pos).at(q0, aa);
          }
        } else {
          if (q1 == aa) {
            de.at(pos, aa) = 0;
          } else if (q0 == aa) {
            de.at(pos, aa) = -tmp;
          } else {
            de.at(pos, aa) += params->h.at(q1, pos) - params->h.at(q0, pos);
            for (size_t pos2 = 0; pos2 < n; pos2++) {
              if (pos2 < i) {
                de.at(pos, aa) +=
                  params->J.at(pos2, i).at(conf.at(pos2), q1) -
                  params->J.at(pos2, i).at(conf.at(pos2), q0);
              } else if (pos2 > i) {
                de.at(pos, aa) +=
                  params->J.at(i, pos2).at(q1, conf.at(pos2)) -
                  params->J.at(i, pos2).at(q0, conf.at(pos2));
              }
            }
          }
        }
      }
    }
  }

  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {

      if (mode == "sqrt") {
        g = arma::exp(de * -1.0 / 2.0 / temperature);
        lambda = arma::accu(g) - n; // n*exp(0) needs to be subtracted.
      } else if (mode == "tanh") {
        g = .5 + .5 * arma::tanh(de * -1.0 / 2.0 / temperature);
        lambda = arma::accu(g) - .5 * n;
      }
      g = g / lambda;

      double rand = uniform(rng);
      double r_sum = 0.0;
      size_t i = 0;
      size_t q0 = 0;
      size_t q1 = 0;
      bool done = false;
      for (i = 0; i < n; i++) {
        for (q1 = 0; q1 < q; q1++) {
          if (conf.at(i) == q1) {
            continue;
          } else {
            r_sum += g.at(i, q1);
          }
          if (r_sum > rand) {
            done = true;
            break;
          }
        }
        if (done) {
          break;
        }
      }

      double tmp = de.at(i, q1);
      q0 = conf.at(i);
      conf.at(i) = q1;
      for (size_t pos = 0; pos < n; pos++) {
        for (size_t aa = 0; aa < q; aa++) {
          if (pos != i) {
            if (pos < i) {
              de.at(pos, aa) += params->J.at(pos, i).at(conf.at(pos), q1) -
                                params->J.at(pos, i).at(conf.at(pos), q0) -
                                params->J.at(pos, i).at(aa, q1) +
                                params->J.at(pos, i).at(aa, q0);
            } else {
              de.at(pos, aa) += params->J.at(i, pos).at(q1, conf.at(pos)) -
                                params->J.at(i, pos).at(q0, conf.at(pos)) -
                                params->J.at(i, pos).at(q1, aa) +
                                params->J.at(i, pos).at(q0, aa);
            }
          } else {
            if (q1 == aa) {
              de.at(pos, aa) = 0;
            } else if (q0 == aa) {
              de.at(pos, aa) = -tmp;
            } else {
              de.at(pos, aa) += params->h.at(q1, pos) - params->h.at(q0, pos);
              for (size_t pos2 = 0; pos2 < n; pos2++) {
                if (pos2 < i) {
                  de.at(pos, aa) +=
                    params->J.at(pos2, i).at(conf.at(pos2), q1) -
                    params->J.at(pos2, i).at(conf.at(pos2), q0);
                } else if (pos2 > i) {
                  de.at(pos, aa) +=
                    params->J.at(i, pos2).at(q1, conf.at(pos2)) -
                    params->J.at(i, pos2).at(q0, conf.at(pos2));
                }
              }
            }
          }
        }
      }
    }
    for (size_t i = 0; i < n; ++i) {
      (*ptr).at(s, i) = (int)conf.at(i);
    }
  }
  return;
};

ostream&
Graph::print_parameters(ostream& os)
{
  log_out << "printing parameters J" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      for (size_t yi = 0; yi < q; ++yi) {
        for (size_t yj = 0; yj < q; ++yj) {
          os << "J " << i << " " << j << " " << yi << " " << yj << " "
             << params->J.at(i, j).at(yi, yj) << endl;
        }
      }
    }
  }
  log_out << "printing parameters h" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t yi = 0; yi < q; ++yi) {
      os << "h " << i << " " << yi << " " << params->h.at(yi, i) << endl;
    }
  }
  log_out << "done" << endl;
  return os;
};

void
Graph::print_parameters(FILE* of)
{
  log_out << "printing parameters J" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      for (size_t yi = 0; yi < q; ++yi) {
        for (size_t yj = 0; yj < q; ++yj) {
          fprintf(of,
                  "J %lu %lu %lu %lu %g\n",
                  i,
                  j,
                  yi,
                  yj,
                  params->J.at(i, j).at(yi, yj));
        }
      }
    }
  }
  log_out << "printing parameters h" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t yi = 0; yi < q; ++yi) {
      fprintf(of, "h %lu %lu %g\n", i, yi, params->h.at(yi, i));
    }
  }
  log_out << "done" << endl;
};
