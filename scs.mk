UNAME = $(shell uname -s)
CC = gcc

ifeq ($(UNAME), Linux)
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS = -lm -lrt
else
# we're on apple, no need to link rt library
LDFLAGS = -lm 
endif

CFLAGS = -g -Wall -pedantic -O3 -funroll-loops -Wstrict-prototypes -Iinclude

LINSYS = linsys
DIRSRC = $(LINSYS)/direct
DIRSRCEXT = $(DIRSRC)/external
INDIRSRC = $(LINSYS)/indirect

OUT = out
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

########### OPTIONAL FLAGS ##########
# CFLAGS += -DDLONG # use longs rather than ints
# CFLAGS += -DFLOAT # use floats rather than doubles
# CFLAGS += -DNOVALIDATE # remove data validation step
# CFLAGS += -DEXTRAVERBOSE # extra verbosity level

############ OPENMP: ############
# uncomment below to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 4 by entering the command:
# export OMP_NUM_THREADS=4

# USE_OPENMP = 1

ifdef USE_OPENMP
    CFLAGS += -fopenmp -DOPENMP
endif

############ SDPS: BLAS + LAPACK ############
# uncomment the line below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

# USE_LAPACK = 1

ifdef USE_LAPACK
  # edit these for your setup:
  LDFLAGS += -lblas -llapack -lgfortran
  CFLAGS += -DLAPACK_LIB_FOUND
# CFLAGS += -DBLAS64 # if blas/lapack lib uses long rather than int
endif
