bin_PROGRAMS = bmdca bmdca_sample

CXX = g++
CXXFLAGS = -O3 -std=c++11 $(OPENMP_CXXFLAGS) $(ARMADILLO_CFLAGS)
LDFLAGS = -lm $(ARMADILLO_LIBS)

DISTCLEANFILES = Makefile.in

bmdca_SOURCES = bmdca.cpp \
                model.cpp \
                msa.cpp \
                msa_stats.cpp \
                run.cpp \
                mcmc.cpp \
                mcmc_stats.cpp \
                graph.cpp \
                utils.cpp

bmdca_sample_SOURCES = bmdca_sample.cpp \
                       generator.cpp \
                       mcmc.cpp \
                       mcmc_stats.cpp \
                       graph.cpp \
                       utils.cpp
