#include "generator.hpp"

Generator::Generator(potts_model params, int n, int q, std::string config_file)
  : N(n)
  , Q(q)
  , graph(n, q)
{
  if (config_file.length() == 0) {
    initializeParameters();
  } else {
    loadParameters(config_file);
  }
  graph.load(params);
};

void
Generator::initializeParameters(void) {
  random_seed = 1;
  adapt_up = 1.5;
  adapt_down = 0.6;
  check_ergo = true;
  t_wait = 10000;
  delta_t = 100;
  temperature = 1.0;
  use_indep_samples = true;
  output_numerical = false;
};

void
Generator::loadParameters(std::string file_name) {
  std::ifstream file(file_name);
  bool reading_bmdca_section = false;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
        continue;
      } else if (line[0] == '[') {
        if (line == "[sampling]") {
          reading_bmdca_section = true;
          continue;
        } else {
          reading_bmdca_section = false;
          continue;
        }
      }
      if (reading_bmdca_section) {
        auto delim_pos = line.find("=");
        auto key = line.substr(0, delim_pos);
        auto value = line.substr(delim_pos + 1);
        setParameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

void
Generator::setParameter(std::string key, std::string value) {
  if (key == "random_seed") {
    random_seed = std::stoi(value);
  } else if (key == "adapt_up") {
    adapt_up = std::stod(value);
  } else if (key == "adapt_down") {
    adapt_down = std::stod(value);
  } else if (key == "check_ergo") {
    if (value.size() == 1) {
      check_ergo = (std::stoi(value) == 1);
    } else {
      check_ergo = (value == "true");
    }
  } else if (key == "t_wait") {
    t_wait = std::stoi(value);
  } else if (key == "delta_t") {
    delta_t = std::stoi(value);
  // } else if (key == "M") {
  //   M = std::stoi(value);
  } else if (key == "temperature") {
    temperature = std::stod(value);
  } else if (key == "use_indep_samples") {
    if (value.size() == 1) {
      use_indep_samples = (std::stoi(value) == 1);
    } else {
      use_indep_samples = (value == "true");
    }
  } else if (key == "output_numerical") {
    if (value.size() == 1) {
      output_numerical = (std::stoi(value) == 1);
    } else {
      output_numerical = (value == "true");
    }
  }
};

void
Generator::sample(arma::Cube<int>* ptr) {
  int M = ptr->n_rows;
  int N = ptr->n_cols;
  int reps = ptr->n_slices;

  if (use_indep_samples) {
#pragma omp parallel
    {
#pragma omp for
      for (int rep = 0; rep < reps; rep++) {
        graph.sample_mcmc((arma::Mat<int>*)&(*ptr).slice(rep),
                          M,
                          t_wait,
                          delta_t,
                          random_seed + rep,
                          temperature);
      }
    }
  }
};

