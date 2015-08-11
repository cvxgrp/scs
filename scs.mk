UNAME = $(shell uname -s)
CC = gcc

ifneq (, $(findstring CYGWIN, $(UNAME)))
ISWINDOWS := 1
else ifneq (, $(findstring MINGW, $(UNAME)))
ISWINDOWS := 1
else ifneq (, $(findstring MSYS, $(UNAME)))
ISWINDOWS := 1
else
ISWINDOWS := 0
endif

ifeq ($(UNAME), Darwin)
# we're on apple, no need to link rt library
LDFLAGS += -lm
SHARED = dylib
else ifeq ($(ISWINDOWS), 1)
# we're on windows (cygwin or msys)
LDFLAGS += -lm
SHARED = dll
else
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS += -lm -lrt
SHARED = so
endif

# Add on default CFLAGS
CFLAGS += -g -Wall -pedantic -O3 -funroll-loops -Wstrict-prototypes -I. -Iinclude
ifneq ($(ISWINDOWS), 1)
CFLAGS += -fPIC
endif

LINSYS = linsys
DIRSRC = $(LINSYS)/direct
INDIRSRC = $(LINSYS)/indirect

OUT = out
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

########### OPTIONAL FLAGS ##########
# these can all be override from the command line
# e.g. make DLONG=1 will override the setting below
DLONG = 0
ifneq ($(DLONG), 0)
CFLAGS += -DDLONG=$(DLONG) # use longs rather than ints
endif
CTRLC = 1
ifneq ($(CTRLC), 0)
CFLAGS += -DCTRLC=$(CTRLC) # graceful interrupts with ctrl-c
endif
FLOAT = 0
ifneq ($(FLOAT), 0)
CFLAGS += -DFLOAT=$(FLOAT) # use floats rather than doubles
endif
NOVALIDATE = 0
ifneq ($(NOVALIDATE), 0)
CFLAGS += -DNOVALIDATE=$(NOVALIDATE)$ # remove data validation step
endif
NOTIMER = 0
ifneq ($(NOTIMER), 0)
CFLAGS += -DNOTIMER=$(NOTIMER) # no timing, times reported as nan
endif
COPYAMATRIX = 1
ifneq ($(COPYAMATRIX), 0)
CFLAGS += -DCOPYAMATRIX=$(COPYAMATRIX) # if normalize, copy A
endif

### VERBOSITY LEVELS: 0,1,2
EXTRAVERBOSE = 0
ifneq ($(EXTRAVERBOSE), 0)
CFLAGS += -DEXTRAVERBOSE=$(EXTRAVERBOSE) # extra verbosity level
endif

############ OPENMP: ############
# set USE_OPENMP = 1 to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 4 by entering the command:
# export OMP_NUM_THREADS=4

USE_OPENMP = 0
ifneq ($(USE_OPENMP), 0)
  CFLAGS += -fopenmp -DOPENMP
  LDFLAGS += -lgomp
endif

############ SDPS: BLAS + LAPACK ############
# set USE_LAPACK = 1 below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

USE_LAPACK = 1
ifneq ($(USE_LAPACK), 0)
  # edit these for your setup:
  BLASLDFLAGS = -lblas -llapack #-lgfortran
  LDFLAGS += $(BLASLDFLAGS)
  CFLAGS += -DLAPACK_LIB_FOUND

  BLAS64 = 0
  ifneq ($(BLAS64), 0)
  CFLAGS += -DBLAS64=$(BLAS64) # if blas/lapack lib uses 64 bit ints
  endif

  NOBLASSUFFIX = 0
  ifneq ($(NOBLASSUFFIX), 0)
  CFLAGS += -DNOBLASSUFFIX=$(NOBLASSUFFIX) # hack to strip blas suffix
  endif

  BLASSUFFIX = "_"
  ifneq ($(BLASSUFFIX), "_")
  CFLAGS += -DBLASSUFFIX=$(BLASSUFFIX) # blas suffix (underscore usually)
  endif
endif

MATLAB_MEX_FILE = 0
ifneq ($(MATLAB_MEX_FILE), 0)
CFLAGS += -DMATLAB_MEX_FILE=$(MATLAB_MEX_FILE) # matlab mex
endif
PYTHON = 0
ifneq ($(PYTHON), 0)
CFLAGS += -DPYTHON=$(PYTHON) # python extension
endif
USING_R = 0
ifneq ($(USING_R), 0)
CFLAGS += -DUSING_R=$(USING_R) # R extension
endif

# debug to see var values, e.g. 'make print-OBJECTS' shows OBJECTS value
print-%: ; @echo $*=$($*)
