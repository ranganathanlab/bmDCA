bin_PROGRAMS = bmdca

CXX = g++
CXXFLAGS = -O3 -std=c++11 -fopenmp $(ARMADILLO_CFLAGS)
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
