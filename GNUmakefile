###############################################################################
#                               GNUmakefile                                   #
#                                                                             #
#   The targets of the whole repository, read by GNU make before the file     #
#   makefile, which stays the makefile of the library that the makefiles of   #
#   examples/, test/ and tools/ include.                                      #
#                                                                             #
#   make             builds the examples, the tools and the tester            #
#   make test        runs the cases of test/tests.txt (MPI apart)             #
#   make test_mpi    runs the cases with MPI under mpiexec                    #
#   make bench       builds the benchmark tools                               #
#   make install     builds lib/libPBeFS.a and copies it, with the            #
#                    headers, into $(PREFIX)/lib and $(PREFIX)/include        #
#                    [see makefile-lib]                                       #
#                                                                             #
#   CC chooses the compiler (g++ by default, e.g., CC=clang++), for the       #
#   library as for the programs.                                              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
###############################################################################

PREFIX ?= /usr/local

all:
	$(MAKE) -C examples
	$(MAKE) -C tools
	$(MAKE) -C test

test:
	$(MAKE) -C test check

test_mpi:
	$(MAKE) -C test check_mpi

bench:
	$(MAKE) -C tools release mpi

# the library [see makefile-lib]

lib:
	$(MAKE) -f makefile-lib

install:
	$(MAKE) -f makefile-lib install PREFIX=$(PREFIX)

clean::
	$(MAKE) -C examples clean
	$(MAKE) -C tools clean
	$(MAKE) -C test clean
	$(MAKE) -f makefile-lib clean
	rm -f obj/*.o

.PHONY: all test test_mpi bench lib install clean

############################# End of GNUmakefile ##############################
