ifeq ($(OS),Windows_NT)
UNAME = CYGWINorMINGWorMSYS
else
UNAME = $(shell uname -s)
endif

#CC = gcc
# For cross-compiling with mingw use these.
#CC = i686-w64-mingw32-gcc -m32
#CC = x86_64-w64-mingw32-gcc-4.8
CUCC = $(CC) #Don't need to use nvcc, since using cuda blas APIs

# For GPU must add cuda libs to path, e.g.
# export DYLD_LIBRARY_PATH=/usr/local/cuda/lib:$DYLD_LIBRARY_PATH

ifneq (, $(findstring CYGWIN, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring MINGW, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring MSYS, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring mingw, $(CC)))
ISWINDOWS := 1
else
ISWINDOWS := 0
endif
endif
endif
endif

ifeq ($(UNAME), Darwin)
# we're on apple, no need to link rt library
LDFLAGS += -lm
SHARED = dylib
SONAME = -install_name
else
ifeq ($(ISWINDOWS), 1)
# we're on windows (cygwin or msys)
LDFLAGS += -lm
SHARED = dll
SONAME = -soname
else
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS += -lm -lrt
SHARED = so
SONAME = -soname
endif
endif

#TODO: check if this works for all platforms:
ifeq ($(CUDA_PATH), )
CUDA_PATH=/usr/local/cuda
endif
CULDFLAGS = -L$(CUDA_PATH)/lib -L$(CUDA_PATH)/lib64 -lcudart -lcublas -lcusparse
CUDAFLAGS = $(CFLAGS) -I$(CUDA_PATH)/include -Ilinsys/gpu -Wno-c++11-long-long # turn off annoying long-long warnings in cuda header files

# Add on default CFLAGS
OPT = -O3
override CFLAGS += -g -Wall -Wwrite-strings -pedantic -funroll-loops -Wstrict-prototypes -I. -Iinclude -Ilinsys $(OPT)
ifneq ($(ISWINDOWS), 1)
override CFLAGS += -fPIC
endif

LINSYS = linsys
DIRSRC = $(LINSYS)/cpu/direct
INDIRSRC = $(LINSYS)/cpu/indirect
GPUDIR = $(LINSYS)/gpu/direct
GPUINDIR = $(LINSYS)/gpu/indirect
MKLSRC = $(LINSYS)/mkl/direct

EXTSRC = $(LINSYS)/external

OUT = out
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib
INSTALL = install

ifeq ($(PREFIX),)
  PREFIX = /usr/local
endif

########### CUSTOM FLAGS ##########
# these can all be override from the command line
CUSTOM_FLAGS =
# e.g. make DLONG=1 will override the setting below
DLONG = 0
ifneq ($(DLONG), 0)
CUSTOM_FLAGS += -DDLONG=$(DLONG) # use longs rather than ints
endif
CTRLC = 1
ifneq ($(CTRLC), 0)
CUSTOM_FLAGS += -DCTRLC=$(CTRLC) # graceful interrupts with ctrl-c
endif
SFLOAT = 0
ifneq ($(SFLOAT), 0)
CUSTOM_FLAGS += -DSFLOAT=$(SFLOAT) # use floats rather than doubles
endif
GPU_TRANSPOSE_MAT = 1
ifneq ($(GPU_TRANSPOSE_MAT), 0)
CUSTOM_FLAGS += -DGPU_TRANSPOSE_MAT=$(GPU_TRANSPOSE_MAT) # transpose A mat in GPU memory
endif
NO_TIMER = 0
ifneq ($(NO_TIMER), 0)
CUSTOM_FLAGS += -DNO_TIMER=$(NO_TIMER) # no timing, times reported as nan
endif
NO_VALIDATE = 0
ifneq ($(NO_VALIDATE), 0)
CUSTOM_FLAGS += -DNO_VALIDATE=$(NO_VALIDATE) # perform problem validation or skip
endif
NO_PRINTING = 0
ifneq ($(NO_PRINTING), 0)
CUSTOM_FLAGS += -DNO_PRINTING=$(NO_PRINTING) # disable printing
endif
NO_READ_WRITE = 0
ifneq ($(NO_READ_WRITE), 0)
CUSTOM_FLAGS += -DNO_READ_WRITE=$(NO_READ_WRITE) # disable printing
endif
### VERBOSITY LEVELS: 0,1,2,...
VERBOSITY = 0
ifneq ($(VERBOSITY), 0)
CUSTOM_FLAGS += -DVERBOSITY=$(VERBOSITY) # verbosity level
endif
COVERAGE = 0
ifneq ($(COVERAGE), 0)
CUSTOM_FLAGS += --coverage # generate test coverage data
endif

# See: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
# This is probably not correct for other systems. TODO: update this
# to work for all combinations of platform / compiler / threading options.
MKLFLAGS = -L$(MKLROOT) -L$(MKLROOT)/lib -Wl,--no-as-needed -lmkl_rt -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -ldl

############ OPENMP: ############
# set USE_OPENMP = 1 to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 4 by entering the command:
# export OMP_NUM_THREADS=4

USE_OPENMP = 0
ifneq ($(USE_OPENMP), 0)
  override CFLAGS += -fopenmp
  LDFLAGS += -fopenmp
endif

############ SDPS: BLAS + LAPACK ############
# set USE_LAPACK = 1 below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

BLASLDFLAGS =
USE_LAPACK = 1
ifneq ($(USE_LAPACK), 0)
  # edit these for your setup:
  BLASLDFLAGS += -llapack -lblas # -lgfortran
  CUSTOM_FLAGS += -DUSE_LAPACK

  BLAS64 = 0
  ifneq ($(BLAS64), 0)
  CUSTOM_FLAGS += -DBLAS64=$(BLAS64) # if blas/lapack lib uses 64 bit ints
  endif

  NOBLASSUFFIX = 0
  ifneq ($(NOBLASSUFFIX), 0)
  CUSTOM_FLAGS += -DNOBLASSUFFIX=$(NOBLASSUFFIX) # hack to strip blas suffix
  endif

  BLASSUFFIX = "_"
  ifneq ($(BLASSUFFIX), "_")
  CUSTOM_FLAGS += -DBLASSUFFIX=$(BLASSUFFIX) # blas suffix (underscore usually)
  endif
endif

MATLAB_MEX_FILE = 0
ifneq ($(MATLAB_MEX_FILE), 0)
CUSTOM_FLAGS += -DMATLAB_MEX_FILE=$(MATLAB_MEX_FILE) # matlab mex
endif
PYTHON = 0
ifneq ($(PYTHON), 0)
CUSTOM_FLAGS += -DPYTHON=$(PYTHON) # python extension
endif
USING_R = 0
ifneq ($(USING_R), 0)
CUSTOM_FLAGS += -DUSING_R=$(USING_R) # R extension
endif

# debug to see var values, e.g. 'make print-OBJECTS' shows OBJECTS value
print-%: ; @echo $*=$($*)

override CFLAGS += $(CUSTOM_FLAGS)
CUDAFLAGS += $(CUSTOM_FLAGS)
