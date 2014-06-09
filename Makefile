# MAKEFILE for scs
include scs.mk

OBJECTS = src/scs.o src/util.o src/cones.o src/cs.o src/linAlg.o

SRC_FILES = $(wildcard src/*.c)
INC_FILES = $(wildcard include/*.h)

AMD_SOURCE = $(wildcard $(DIRSRCEXT)/amd_*.c)
DIRECT_OBJECTS = $(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = $(OUT)/demo_direct $(OUT)/demo_indirect $(OUT)/demo_SOCP_indirect $(OUT)/demo_SOCP_direct

.PHONY: default 

default: $(TARGETS) $(OUT)/libscsdir.a $(OUT)/libscsindir.a 
    #$(OUT)/libscsdir.so $(OUT)/libscsindir.so 	
	@echo "**********************************************************************************"
	@echo "Successfully compiled scs, copyright Brendan O'Donoghue 2014."
	@echo "To test, type '$(OUT)/demo_direct' or '$(OUT)/demo_indirect'."
	@echo "**********************************************************************************"
ifneq ($(USE_LAPACK), 0)
	@echo "Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, and EXPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, and EXPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to point to the library"
	@echo "install locations, and recompile with 'make purge', 'make'."
endif
	@echo "**********************************************************************************"

src/scs.o	: $(SRC_FILES) $(INC_FILES)
src/util.o	: src/util.c include/util.h
src/cones.o	: src/cones.c include/cones.h
src/cs.o	: src/cs.c include/cs.h
src/linAlg.o: src/linAlg.c include/linAlg.h

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

#$(OUT)/libscsdir.so: $(OBJECTS) $(DIRSRC)/private.o $(DIRECT_OBJECTS) $(LINSYS)/common.o
#	mkdir -p $(OUT)
#	$(CC) -shared -o $@ $^ $(LDFLAGS)

#$(OUT)/libscsindir.so: $(OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/common.o
#	mkdir -p $(OUT)
#	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(OUT)/demo_direct: examples/c/demo.c $(OUT)/libscsdir.a
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" $^ -o $@ $(LDFLAGS)

$(OUT)/demo_indirect: examples/c/demo.c $(OUT)/libscsindir.a
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" $^  -o $@ $(LDFLAGS)

$(OUT)/demo_SOCP_direct: examples/c/randomSOCPProb.c $(OUT)/libscsdir.a
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/demo_SOCP_indirect: examples/c/randomSOCPProb.c $(OUT)/libscsindir.a
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

## To compile dense versions: make dense
## Need to connect to blas lib (set USE_LAPACK = 1 in scs.mk)

DENSE=linsys_dense
DENSE_OUT=$(DENSE)/out
DENSE_TARGETS=$(DENSE)/direct/private.o $(DENSE)/indirect/private.o $(DENSE)/common.o $(DENSE)/libscsdir.a $(DENSE)/libscsindir.a $(DENSE_OUT)/demo_SOCP_direct $(DENSE_OUT)/demo_SOCP_indirect

dense: $(DENSE_TARGETS)

$(DENSE)/direct/private.o: $(DENSE)/direct/private.c  $(DENSE)/direct/private.h
$(DENSE)/indirect/private.o: $(DENSE)/indirect/private.c $(DENSE)/indirect/private.h
$(DENSE)/common.o: $(DENSE)/common.c $(DENSE)/common.h

$(DENSE)/libscsdir.a: $(OBJECTS) $(DENSE)/direct/private.o $(DENSE)/common.o
	$(ARCHIVE) $(DENSE)/libscsdir.a $^
	- $(RANLIB) $(DENSE)/libscsdir.a

$(DENSE)/libscsindir.a: $(OBJECTS) $(DENSE)/indirect/private.o $(DENSE)/common.o
	$(ARCHIVE) $(DENSE)/libscsindir.a $^
	- $(RANLIB) $(DENSE)/libscsindir.a

$(DENSE_OUT)/demo_SOCP_direct: $(DENSE)/randomSOCPProb.c $(DENSE)/libscsdir.a
	mkdir -p $(DENSE_OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(DENSE_OUT)/demo_SOCP_indirect: $(DENSE)/randomSOCPProb.c $(DENSE)/libscsindir.a
	mkdir -p $(DENSE_OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) $(LINSYS)/common.o $(DIRSRC)/private.o $(INDIRSRC)/private.o $(DENSE_TARGETS)
	@rm -rf $(OUT)/*.dSYM
	@rm -rf matlab/*.mex*
	@rm -rf .idea
	@rm -rf python/*.pyc
	@rm -rf python/build
purge: clean 
	@rm -rf $(OUT)

