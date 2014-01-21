UNAME = $(shell uname -s)
CC = gcc
CFLAGS = -g -Wall -pedantic -O3 -Isrc/include -funroll-loops
LDFLAGS = -lm 
# for openmp:
CFLAGS += -fopenmp

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

ifeq (1,1)
  CFLAGS += -I/opt/local/include -I/usr/local/include
  LDFLAGS += -L/opt/local/lib -L/usr/local/lib
  LDFLAGS += -lopenblas -llapack -llapacke
  CFLAGS += -DLAPACK_LIB_FOUND
endif

AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib
