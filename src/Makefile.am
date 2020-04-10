bin_PROGRAMS = bmdca bmdca_sample arma2ascii

CXXFLAGS = -O3 -std=c++11 $(OPENMP_CXXFLAGS) $(ARMADILLO_CFLAGS) \
           -DARMA_DONT_USE_HDF5
if LINUX
AM_LDFLAGS = -lm $(ARMADILLO_LIBS) -lopenblas -llapack
endif
if MACOS
AM_LDFLAGS = -lm $(ARMADILLO_LIBS) -framework Accelerate
endif
if WINDOWS
AM_LDFLAGS = -static -static-libgcc -static-libstdc++ \
             -lm $(ARMADILLO_LIBS).dll -lopenblas -llapack
endif

DISTCLEANFILES = Makefile.in

bmdca_SOURCES = bmdca.cpp \
                model.cpp \
                msa.cpp \
                msa_stats.cpp \
                run.cpp \
                mcmc.cpp \
                mcmc_stats.cpp \
                graph_arma.cpp \
                utils.cpp

bmdca_sample_SOURCES = bmdca_sample.cpp \
                       generator.cpp \
                       mcmc.cpp \
                       mcmc_stats.cpp \
                       graph_arma.cpp \
                       utils.cpp

arma2ascii_SOURCES = arma_convert.cpp \
                     utils.cpp
