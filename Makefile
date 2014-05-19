# MAKEFILE for scs
include scs.mk

OBJECTS = src/scs.o src/util.o src/cones.o src/cs.o src/linAlg.o
AMD_SOURCE = $(wildcard $(DIRSRCEXT)/amd_*.c)
DIRECT_OBJECTS = $(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = $(OUT)/demo_direct $(OUT)/demo_indirect $(OUT)/demo_SOCP_indirect $(OUT)/demo_SOCP_direct

.PHONY: default 
default: $(OUT)/libscsdir.a $(OUT)/libscsindir.a $(TARGETS)
	@echo "**********************************************************************************"
	@echo "Successfully compiled scs, copyright Brendan O'Donoghue 2014."
	@echo "To test, type '$(OUT)/demo_direct' or '$(OUT)/demo_indirect'."
	@echo "**********************************************************************************"
ifdef USE_LAPACK
	@echo "Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, and EXPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, and EXPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to point to the library"
	@echo "install locations, and recompile with 'make purge', 'make'."
endif
	@echo "**********************************************************************************"

scs 	: src/scs.c include/scs.h include/linSys.h include/normalize.h include/util.h
util	: src/util.c include/util.h
cones	: src/cones.c include/cones.h
cs		: src/cs.c include/cs.h
linAlg  : src/linAlg.c include/linAlg.h

$(DIRSRCEXT)/private.o		    : $(DIRSRC)/private.h $(LINSYS)/common.h
$(DIRSRCEXT)/ldl.o			    : $(DIRSRCEXT)/ldl.h
$(DIRSRCEXT)/amd_1.o			: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_2.o			: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_aat.o		    : $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_control.o	    : $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_defaults.o 	: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_dump.o		    : $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_global.o		: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_info.o		    : $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_order.o		: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_post_tree.o	: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_postorder.o	: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_preprocess.o	: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(DIRSRCEXT)/amd_valid.o		: $(DIRSRCEXT)/amd_internal.h $(DIRSRCEXT)/amd.h
$(INDIRSRC)/private.o	        : $(INDIRSRC)/private.h $(LINSYS)/common.h
$(LINSYS)/common.o              : $(LINSYS)/common.h

$(OUT)/libscsdir.a: $(OBJECTS) $(DIRSRC)/private.o $(DIRECT_OBJECTS) $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(ARCHIVE) $(OUT)/libscsdir.a $^
	- $(RANLIB) $(OUT)/libscsdir.a

$(OUT)/libscsindir.a: $(OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/common.o
	mkdir -p $(OUT)
	$(ARCHIVE) $(OUT)/libscsindir.a $^
	- $(RANLIB) $(OUT)/libscsindir.a

$(OUT)/demo_direct: examples/c/demo.c $(OUT)/libscsdir.a examples/c/problemUtils.h
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" -o $@ $^ $(LDFLAGS) 

$(OUT)/demo_indirect: examples/c/demo.c $(OUT)/libscsindir.a examples/c/problemUtils.h
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" -o $@ $^ $(LDFLAGS) 

$(OUT)/demo_SOCP_direct: examples/c/randomSOCPProb.c $(OUT)/libscsdir.a examples/c/problemUtils.h
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/demo_SOCP_indirect: examples/c/randomSOCPProb.c $(OUT)/libscsindir.a examples/c/problemUtils.h
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

## To compile dense versions: make dense
## Need to connect to blas lib (set USE_LAPACK = 1 in scs.mk)

DENSE=linsys_dense
DENSE_OUT=$(DENSE)/out
DENSE_TARGETS=$(DENSE)/direct/private.o $(DENSE)/indirect/private.o $(DENSE)/common.o $(DENSE)/libscsdir.a $(DENSE)/libscsindir.a $(DENSE_OUT)/demo_SOCP_direct $(DENSE_OUT)/demo_SOCP_indirect

dense: $(DENSE_TARGETS)

$(DENSE)/direct/private.o: $(DENSE)/direct/private.h
$(DENSE)/indirect/private.o: $(DENSE)/indirect/private.h
$(DENSE)/common.o: $(DENSE)/common.h

$(DENSE)/libscsdir.a: $(OBJECTS) $(DENSE)/direct/private.o $(DENSE)/common.o
	$(ARCHIVE) $(DENSE)/libscsdir.a $^
	- $(RANLIB) $(DENSE)/libscsdir.a

$(DENSE)/libscsindir.a: $(OBJECTS) $(DENSE)/indirect/private.o $(DENSE)/common.o
	$(ARCHIVE) $(DENSE)/libscsindir.a $^
	- $(RANLIB) $(DENSE)/libscsindir.a

$(DENSE_OUT)/demo_SOCP_direct: $(DENSE)/randomSOCPProb.c $(DENSE)/libscsdir.a $(DENSE)/problemUtils.h
	mkdir -p $(DENSE_OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(DENSE_OUT)/demo_SOCP_indirect: $(DENSE)/randomSOCPProb.c $(DENSE)/libscsindir.a $(DENSE)/problemUtils.h
	mkdir -p $(DENSE_OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) $(DIRSRC)/private.o $(INDIRSRC)/private.o $(DENSE_TARGETS)
	@rm -rf $(OUT)/*.dSYM
	@rm -rf matlab/*.mex*
	@rm -rf .idea
	@rm -rf python/*.pyc
	@rm -rf python/build
purge: clean 
	@rm -rf $(OUT)

