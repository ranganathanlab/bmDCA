## BEGIN HEADER ##
tt1 = gasdev
tt2 = graph1
tt3 = graph2
tt4 = graph3
tt5 = graphs
tt6 = do_montecarlo
## END HEADER ##

## BEGIN OTHERTAGS ##
## END OTHERTAGS ##

## BEGIN COMPILERS ##
COMPILER_C = gcc
COMPILER_CPP = g++
LINKER_C = gcc
LINKER_CPP = g++
## END COMPILERS ##

## BEGIN FLAGS ##
FLAGS = -Wall -g -O3 -fopenmp -I..
FLAGSL = -lm -fopenmp
## END FLAGS ##

## BEGIN OBJECTS ##
tt1obj = $(tt1).o
tt2obj = $(tt2).o
tt3obj = $(tt3).o
tt4obj = $(tt4).o
tt5obj = $(tt5).o
tt6obj = $(tt6).o $(tt5).o $(tt4).o $(tt3).o $(tt2).o $(tt1).o
## END OBJECTS ##

## BEGIN OUTPUTS ##
tt6out = ./do_montecarlo
## END OUTPUTS ##

## BEGIN DEPENDENCIES ##
tt1deps = $(tt1).cpp $(tt1).hpp
tt2deps = $(tt2).cpp $(tt2).hpp $(tt1).hpp $(tt5).hpp $(tt3).hpp $(tt4).hpp
tt3deps = $(tt3).cpp $(tt3).hpp $(tt1).hpp $(tt2).hpp $(tt4).hpp $(tt5).hpp
tt4deps = $(tt4).cpp $(tt4).hpp $(tt1).hpp $(tt2).hpp $(tt3).hpp $(tt5).hpp
tt5deps = $(tt5).cpp $(tt5).hpp
tt6deps = $(tt6).cpp $(tt5).hpp $(tt2).hpp $(tt3).hpp $(tt4).hpp
## END DEPENDENCIES ##

## BEGIN PHONYS ##
.PHONY : all clean ctags tt6

all : $(tt6out)
## END PHONYS ##

## BEGIN BODY tt1 ##
$(tt1).o : $(tt1deps)
	$(COMPILER_CPP) -c $(FLAGS) $(tt1).cpp
## END BODY tt1 ##

## BEGIN BODY tt2 ##
$(tt2).o : $(tt2deps)
	$(COMPILER_CPP) -c $(FLAGS) $(tt2).cpp
## END BODY tt2 ##

## BEGIN BODY tt3 ##
$(tt3).o : $(tt3deps)
	$(COMPILER_CPP) -c $(FLAGS) $(tt3).cpp
## END BODY tt3 ##

## BEGIN BODY tt4 ##
$(tt4).o : $(tt4deps)
	$(COMPILER_CPP) -c $(FLAGS) $(tt4).cpp
## END BODY tt4 ##

## BEGIN BODY tt5 ##
$(tt5).o : $(tt5deps)
	$(COMPILER_CPP) -c $(FLAGS) $(tt5).cpp
## END BODY tt5 ##

## BEGIN BODY tt6 ##
tt6 : $(tt6out)

$(tt6out) : $(tt6obj)
	$(LINKER_CPP) $(FLAGSL) -o $(tt6out) $(tt6obj)

$(tt6).o : $(tt6deps)
	$(COMPILER_CPP) -c $(FLAGS) $(tt6).cpp
## END BODY tt6 ##

## BEGIN CLEAN ##
clean :
	rm -f ctags *~ core.* $(tt1obj) $(tt2obj) $(tt3obj) $(tt4obj) $(tt5obj) $(tt6obj)
## END CLEAN ##

## BEGIN TAGS ##
ctags : tags

tags : $(tt1deps) $(tt2deps) $(tt3deps) $(tt4deps) $(tt5deps) $(tt6deps)
	ctags -R
## END TAGS ##

## BEGIN USERDEFINED ##
## END USERDEFINED ##
