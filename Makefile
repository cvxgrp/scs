# MAKEFILE for scs
include scs.mk

OBJECTS = src/scs.o src/util.o src/cones.o src/cs.o src/linAlg.o
AMD_SOURCE = $(wildcard direct/dirsrc/amd_*.c)
DIRECT_OBJECTS = direct/dirsrc/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = bin/demo_direct bin/demo_indirect

.PHONY: default 
default: lib/libscsdir.a lib/libscsindir.a bin/demo_direct bin/demo_indirect
	@echo "**********************************************************************************"
ifdef USE_LAPACK
	@echo "Compiled with blas and lapack, can solve LPs, SOCPS, SDPs, and EXPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, and EXPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to point to the library"
	@echo "install locations, and recompile with 'make purge', 'make'."
endif
	@echo "**********************************************************************************"

scs 	: src/scs.c src/scs.h src/linAlg.h
util	: src/util.c src/util.h
cones	: src/cones.c src/cones.h
cs		: src/cs.c src/cs.h
linAlg  : src/linAlg.c src/linAlh.h

direct/dirsrc/private.o		: direct/private.h
direct/dirsrc/ldl.o			: direct/dirsrc/ldl.h
direct/dirsrc/amd_1.o			: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_2.o			: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_aat.o		: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_control.o	: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_defaults.o 	: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_dump.o		: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_global.o		: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_info.o		: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_order.o		: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_post_tree.o	: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_postorder.o	: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_preprocess.o	: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
direct/dirsrc/amd_valid.o		: direct/dirsrc/amd_internal.h direct/dirsrc/amd.h
indirect/dirsrc/private.o	: indirect/private.h

lib/libscsdir.a: $(OBJECTS) direct/private.o  $(DIRECT_OBJECTS)
	mkdir -p lib
	$(ARCHIVE) lib/libscsdir.a $^
	- $(RANLIB) lib/libscsdir.a

lib/libscsindir.a: $(OBJECTS) indirect/private.o
	mkdir -p lib
	$(ARCHIVE) lib/libscsindir.a $^
	- $(RANLIB) lib/libscsindir.a

bin/demo_direct: src/run_scs.c lib/libscsdir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_sparse\"" -o $@ $^ $(LDFLAGS) 

bin/demo_indirect: src/run_scs.c lib/libscsindir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_sparse\"" -o $@ $^ $(LDFLAGS) 

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) direct/private.o indirect/private.o 
	@rm -rf bin/*.dSYM
	@rm -rf matlab/*.mex*
purge: clean 
	@rm -rf bin lib
