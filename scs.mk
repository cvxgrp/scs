UNAME = $(shell uname -s)
CC = gcc

ifeq ($(UNAME), Linux)
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS = -lm -lrt
SHARED = so
else
# we're on apple, no need to link rt library
LDFLAGS = -lm
SHARED = dylib
endif

CFLAGS = -g -Wall -pedantic -O3 -funroll-loops -Wstrict-prototypes -fPIC -I. -Iinclude #-Wextra

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
# CFLAGS += -DNOBLASUNDERSCORE # if your blas install does not use underscores in function names 

############ OPENMP: ############
# set USE_OPENMP = 1 to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 4 by entering the command:
# export OMP_NUM_THREADS=4

USE_OPENMP = 0

ifneq ($(USE_OPENMP), 0)
  CFLAGS += -fopenmp -DOPENMP
# LDFLAGS += -lgomp
endif

############ SDPS: BLAS + LAPACK ############
# set USE_LAPACK = 1 below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

USE_LAPACK = 0

ifneq ($(USE_LAPACK), 0)
  # edit these for your setup:
  LDFLAGS += -lblas -llapack #-lgfortran
  CFLAGS += -DLAPACK_LIB_FOUND
  # CFLAGS += -DBLAS64 # if blas/lapack lib uses long rather than int
endif
