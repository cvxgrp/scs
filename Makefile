# MAKEFILE for scs
include scs.mk

SCS_OBJECTS = src/scs.o src/util.o src/cones.o src/cs.o src/linAlg.o src/ctrlc.o src/scs_version.o

SRC_FILES = $(wildcard src/*.c)
INC_FILES = $(wildcard include/*.h)

AMD_SOURCE = $(wildcard $(DIRSRCEXT)/amd_*.c)
DIRECT_SCS_OBJECTS = $(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o)
TARGETS = $(OUT)/demo_direct $(OUT)/demo_indirect $(OUT)/demo_SOCP_indirect $(OUT)/demo_SOCP_direct

.PHONY: default 

default: $(TARGETS) $(OUT)/libscsdir.a $(OUT)/libscsindir.a $(OUT)/libscsdir.$(SHARED) $(OUT)/libscsindir.$(SHARED)
	@echo "****************************************************************************************"
	@echo "Successfully compiled scs, copyright Brendan O'Donoghue 2012."
	@echo "To test, type '$(OUT)/demo_direct' or '$(OUT)/demo_indirect',"
	@echo "or '$(OUT)/demo_SOCP_indirect' to solve a random SOCP."
	@echo "**********************************************************************************"
ifneq ($(USE_LAPACK), 0)
	@echo "Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, ECPs, and PCPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, ECPs, and PCPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to set USE_LAPACK=1"
	@echo "and point to the library install locations, and recompile with 'make purge', 'make'."
endif
	@echo "****************************************************************************************"

src/scs.o	: $(SRC_FILES) $(INC_FILES)
src/util.o	: src/util.c include/util.h include/constants.h
src/cones.o	: src/cones.c include/cones.h include/scs_blas.h
src/cs.o	: src/cs.c include/cs.h
src/linAlg.o: src/linAlg.c include/linAlg.h
src/ctrl.o  : src/ctrl.c include/ctrl.h
src/scs_version.o: src/scs_version.c include/constants.h

$(DIRSRC)/private.o: $(DIRSRC)/private.c  $(DIRSRC)/private.h
$(INDIRSRC)/indirect/private.o: $(INDIRSRC)/private.c $(INDIRSRC)/private.h
$(LINSYS)/common.o: $(LINSYS)/common.c $(LINSYS)/common.h

$(OUT)/libscsdir.a: $(SCS_OBJECTS) $(DIRSRC)/private.o $(DIRECT_SCS_OBJECTS) $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(ARCHIVE) $(OUT)/libscsdir.a $^
	- $(RANLIB) $(OUT)/libscsdir.a

$(OUT)/libscsindir.a: $(SCS_OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(ARCHIVE) $(OUT)/libscsindir.a $^
	- $(RANLIB) $(OUT)/libscsindir.a

$(OUT)/libscsdir.$(SHARED): $(SCS_OBJECTS) $(DIRSRC)/private.o $(DIRECT_SCS_OBJECTS) $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(OUT)/libscsindir.$(SHARED): $(SCS_OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(OUT)/demo_direct: examples/c/demo.c $(OUT)/libscsdir.a
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" $^ -o $@ $(LDFLAGS)

$(OUT)/demo_indirect: examples/c/demo.c $(OUT)/libscsindir.a
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" $^  -o $@ $(LDFLAGS)

$(OUT)/demo_SOCP_direct: examples/c/randomSOCPProb.c $(OUT)/libscsdir.$(SHARED)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/demo_SOCP_indirect: examples/c/randomSOCPProb.c $(OUT)/libscsindir.$(SHARED)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

NVCC = nvcc -g
GPU = gpu
gpu: $(OUT)/demo_SOCP_indirect_gpu 

$(GPU)/private.o: $(GPU)/private.cu
	$(NVCC) -c -I include/ -I/Developer/NVIDIA/CUDA-7.5/samples/common/inc/ -Xcompiler -O3,-fPIC -o $(GPU)/private.o $^

$(GPU)/libscsindirgpu.a: $(SCS_OBJECTS) $(GPU)/private.o
	$(ARCHIVE) $(GPU)/libscsindirgpu.a $^
	- $(RANLIB) $(GPU)/libscsindirgpu.a

$(OUT)/demo_SOCP_indirect_gpu: examples/c/randomSOCPProb.c $(SCS_OBJECTS) $(GPU)/private.o
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -I/Developer/NVIDIA/CUDA-7.5/samples/common/inc/ -L/usr/local/cuda/lib -lcudart -lcublas -lcusparse -lcudadevrt

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(SCS_OBJECTS) $(DIRECT_SCS_OBJECTS) $(LINSYS)/common.o $(DIRSRC)/private.o $(INDIRSRC)/private.o $(GPU)/private.o
	@rm -rf $(OUT)/*.dSYM
	@rm -rf $(OUT)/*gpu*
	@rm -rf matlab/*.mex*
	@rm -rf .idea
	@rm -rf python/*.pyc
	@rm -rf python/build
purge: clean 
	@rm -rf $(OUT)

