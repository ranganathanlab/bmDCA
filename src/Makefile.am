bin_PROGRAMS = bmdca

CXX = g++
CXXFLAGS = -O3 -std=c++11 -fopenmp
LDFLAGS = -larmadillo -lm

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
