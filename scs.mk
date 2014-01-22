UNAME = $(shell uname -s)
CC = gcc
CFLAGS = -g -Wall -pedantic -O3 -Isrc/include -funroll-loops
LDFLAGS = -lm 

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

############ OPENMP: ############
# uncomment below to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 16 using:
# export OMP_NUM_THREADS=16

# CFLAGS += -fopenmp

############ SDPS: BLAS + LAPACK ############
# uncomment the line below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

# USE_LAPACK = 1

ifdef USE_LAPACK
  # change these for your setup:
  CFLAGS += -I/opt/local/include -I/usr/local/include
  LDFLAGS += -L/opt/local/lib -L/usr/local/lib
  LDFLAGS += -lopenblas -llapack -llapacke 
  CFLAGS += -DLAPACK_LIB_FOUND
endif
