UNAME = $(shell uname -s)
CC = gcc
LDFLAGS = -lm 

CFLAGS = -g -Wall -pedantic -O3 -Iinclude -funroll-loops 

DIRSRC = linsys/direct
DIRSRCEXT = $(DIRSRC)/external
INDIRSRC = linsys/indirect

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
# set the number of threads to, for example, 4 using:
# export OMP_NUM_THREADS=4

# CFLAGS += -fopenmp

############ SDPS: BLAS + LAPACK ############
# uncomment the line below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

# USE_LAPACK = 1

ifdef USE_LAPACK
  # edit these for your setup:
  CFLAGS += -I/opt/local/include -I/usr/local/include
  LDFLAGS += -L/opt/local/lib -L/usr/local/lib
  LDFLAGS += -lopenblas -llapack -llapacke 
  CFLAGS += -DLAPACK_LIB_FOUND
  # lapacke + cblas are not ansi c90 compliant:
  ifeq ($(UNAME), Darwin)
      CFLAGS   += -std=c99
  else
      CFLAGS   += -std=gnu99
  endif
endif
