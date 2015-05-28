# MAKEFILE for scs
include scs.mk

OBJECTS = src/scs.o src/util.o src/cones.o src/cs.o src/linAlg.o src/ctrlc.o src/scs_version.o

SRC_FILES = $(wildcard src/*.c)
INC_FILES = $(wildcard include/*.h)

AMD_SOURCE = $(wildcard $(DIRSRCEXT)/amd_*.c)
DIRECT_OBJECTS = $(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = $(OUT)/demo_direct $(OUT)/demo_indirect $(OUT)/demo_SOCP_indirect $(OUT)/demo_SOCP_direct

.PHONY: default 

default: $(TARGETS) $(OUT)/libscsdir.a $(OUT)/libscsindir.a $(OUT)/libscsdir.$(SHARED) $(OUT)/libscsindir.$(SHARED)
	@echo "**********************************************************************************"
	@echo "Successfully compiled scs, copyright Brendan O'Donoghue 2012."
	@echo "To test, type '$(OUT)/demo_direct' or '$(OUT)/demo_indirect',"
	@echo "or '$(OUT)/demo_SOCP_indirect' to solve a random SOCP."
	@echo "**********************************************************************************"
ifneq ($(USE_LAPACK), 0)
	@echo "Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, and ECPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, and ECPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to set USE_LAPACK=1"
	@echo "and point to the library install locations, and recompile with 'make purge', 'make'."
endif
	@echo "**********************************************************************************"

src/scs.o	: $(SRC_FILES) $(INC_FILES)
src/util.o	: src/util.c include/util.h include/constants.h
src/cones.o	: src/cones.c include/cones.h
src/cs.o	: src/cs.c include/cs.h
src/linAlg.o: src/linAlg.c include/linAlg.h
src/ctrl.o  : src/ctrl.c include/ctrl.h
src/scs_version.o: src/scs_version.c include/constants.h

$(DIRSRC)/private.o: $(DIRSRC)/private.c  $(DIRSRC)/private.h
$(INDIRSRC)/indirect/private.o: $(INDIRSRC)/private.c $(INDIRSRC)/private.h
$(LINSYS)/common.o: $(LINSYS)/common.c $(LINSYS)/common.h

$(OUT)/libscsdir.a: $(OBJECTS) $(DIRSRC)/private.o $(DIRECT_OBJECTS) $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(ARCHIVE) $(OUT)/libscsdir.a $^
	- $(RANLIB) $(OUT)/libscsdir.a

$(OUT)/libscsindir.a: $(OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(ARCHIVE) $(OUT)/libscsindir.a $^
	- $(RANLIB) $(OUT)/libscsindir.a

$(OUT)/libscsdir.$(SHARED): $(OBJECTS) $(DIRSRC)/private.o $(DIRECT_OBJECTS) $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(OUT)/libscsindir.$(SHARED): $(OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/common.o
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

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) $(LINSYS)/common.o $(DIRSRC)/private.o $(INDIRSRC)/private.o
	@rm -rf $(OUT)/*.dSYM
	@rm -rf matlab/*.mex*
	@rm -rf .idea
	@rm -rf python/*.pyc
	@rm -rf python/build
purge: clean 
	@rm -rf $(OUT)

