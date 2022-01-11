# MAKEFILE for scs
include scs.mk

SCS_OBJECTS = src/util.o src/cones.o src/aa.o src/rw.o src/linalg.o src/ctrlc.o src/scs_version.o src/normalize.o
SCS_O = src/scs.o
SCS_INDIR_O = src/scs_indir.o

SRC_FILES = $(wildcard src/*.c)
INC_FILES = $(wildcard include/*.h)

AMD_SOURCE = $(wildcard $(EXTSRC)/amd/*.c)
LDL_SOURCE = $(EXTSRC)/qdldl/qdldl.c
AMD_OBJS = $(AMD_SOURCE:.c=.o)
LDL_OBJS = $(LDL_SOURCE:.c=.o)
TARGETS = $(OUT)/demo_socp_indirect $(OUT)/demo_socp_direct $(OUT)/run_from_file_indirect $(OUT)/run_from_file_direct

.PHONY: default

default: $(TARGETS) $(OUT)/libscsdir.a $(OUT)/libscsindir.a $(OUT)/libscsdir.$(SHARED) $(OUT)/libscsindir.$(SHARED)
	@echo "****************************************************************************************"
	@echo "Successfully compiled scs, copyright Brendan O'Donoghue 2012."
	@echo "To test, type '$(OUT)/demo_socp_direct' to solve a random SOCP."
	@echo "**********************************************************************************"
ifneq ($(USE_LAPACK), 0)
	@echo "Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, ECPs, and PCPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, ECPs, and PCPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to set USE_LAPACK=1"
	@echo "and point to the library install locations, and recompile with 'make purge', 'make'."
endif
	@echo "****************************************************************************************"

$(SCS_O): src/scs.c $(INC_FILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(SCS_INDIR_O): src/scs.c $(INC_FILES)
	$(CC) $(CFLAGS) -DINDIRECT=1 -c $< -o $@

%.o : src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

src/util.o	: src/util.c $(INC_FILES)
src/cones.o	: src/cones.c $(INC_FILES)
src/aa.o	: src/aa.c $(INC_FILES)
src/rw.o	: src/rw.c $(INC_FILES)
src/linalg.o: src/linalg.c $(INC_FILES)
src/ctrl.o  : src/ctrl.c $(INC_FILES)
src/scs_version.o: src/scs_version.c $(INC_FILES)

$(DIRSRC)/private.o: $(DIRSRC)/private.c  $(DIRSRC)/private.h
$(INDIRSRC)/indirect/private.o: $(INDIRSRC)/private.c $(INDIRSRC)/private.h
$(LINSYS)/scs_matrix.o: $(LINSYS)/scs_matrix.c $(LINSYS)/scs_matrix.h
$(LINSYS)/csparse.o: $(LINSYS)/csparse.c $(LINSYS)/csparse.h

$(OUT)/libscsdir.a: $(SCS_O) $(SCS_OBJECTS) $(DIRSRC)/private.o $(AMD_OBJS) $(LDL_OBJS) $(LINSYS)/scs_matrix.o $(LINSYS)/csparse.o
	mkdir -p $(OUT)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

$(OUT)/libscsindir.a: $(SCS_INDIR_O) $(SCS_OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/scs_matrix.o $(LINSYS)/csparse.o
	mkdir -p $(OUT)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

$(OUT)/libscsdir.$(SHARED): $(SCS_O) $(SCS_OBJECTS) $(DIRSRC)/private.o $(AMD_OBJS) $(LDL_OBJS) $(LINSYS)/scs_matrix.o $(LINSYS)/csparse.o
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -shared -Wl,$(SONAME),$(@:$(OUT)/%=%) -o $@ $^ $(LDFLAGS)

$(OUT)/libscsindir.$(SHARED): $(SCS_INDIR_O) $(SCS_OBJECTS) $(INDIRSRC)/private.o $(LINSYS)/scs_matrix.o $(LINSYS)/csparse.o
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -shared -Wl,$(SONAME),$(@:$(OUT)/%=%) -o $@ $^ $(LDFLAGS)

$(OUT)/demo_socp_direct: test/random_socp_prob.c $(OUT)/libscsdir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/demo_socp_indirect: test/random_socp_prob.c $(OUT)/libscsindir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/run_from_file_direct: test/run_from_file.c $(OUT)/libscsdir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/run_from_file_indirect: test/run_from_file.c $(OUT)/libscsindir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/run_from_file_gpu_indirect: test/run_from_file.c $(OUT)/libscsgpuindir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CULDFLAGS)

# basic testing
.PHONY: test
test: $(OUT)/run_tests_indirect $(OUT)/run_tests_direct
$(OUT)/run_tests_indirect: test/run_tests.c $(OUT)/libscsindir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -Itest
$(OUT)/run_tests_direct: test/run_tests.c $(OUT)/libscsdir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -Itest

.PHONY: test_gpu
test_gpu: $(OUT)/run_tests_gpu_indirect # $(OUT)/run_tests_gpu_direct

$(OUT)/run_tests_gpu_indirect: test/run_tests.c $(OUT)/libscsgpuindir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CULDFLAGS) -Itest

# $(OUT)/run_tests_gpu_direct: test/run_tests.c $(OUT)/libscsgpudir.a
# 	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CULDFLAGS) -Itest

# REQUIRES GPU AND CUDA INSTALLED
gpu: gpu_indirect # gpu_direct

# gpu_direct: $(OUT)/demo_socp_gpu_direct $(OUT)/libscsgpudir.$(SHARED) $(OUT)/libscsgpudir.a $(OUT)/run_from_file_gpu_direct
gpu_indirect: $(OUT)/demo_socp_gpu_indirect $(OUT)/libscsgpuindir.$(SHARED) $(OUT)/libscsgpuindir.a $(OUT)/run_from_file_gpu_indirect

$(LINSYS)/gpu/gpu.o: $(LINSYS)/gpu/gpu.c
	$(CUCC) -c -o $@ $^ $(CUDAFLAGS)

# $(GPUDIR)/private.o: $(GPUDIR)/private.c
# 	$(CUCC) -c -o $(GPUDIR)/private.o $^ $(CUDAFLAGS)

$(GPUINDIR)/private.o: $(GPUINDIR)/private.c
	$(CUCC) -c -o $@ $^ $(CUDAFLAGS)

# $(OUT)/libscsgpudir.$(SHARED): $(SCS_O) $(SCS_OBJECTS) $(GPUDIR)/private.o $(AMD_OBJS) $(LINSYS)/scs_matrix.o $(LINSYS)/gpu/gpu.o
#	 mkdir -p $(OUT)
# 	$(CC) $(CFLAGS) -shared -Wl,$(SONAME),$(@:$(OUT)/%=%) -o $@ $^ $(LDFLAGS) $(CULDFLAGS)

# $(OUT)/libscsgpudir.a: $(SCS_INDIR_O) $(SCS_OBJECTS) $(GPUDIR)/private.o $(AMD_OBJS) $(LINSYS)/scs_matrix.o $(LINSYS)/gpu/gpu.o
#  	mkdir -p $(OUT)
# 	$(ARCHIVE) $@ $^
# 	- $(RANLIB) $@

$(OUT)/libscsgpuindir.$(SHARED): $(SCS_O) $(SCS_OBJECTS) $(GPUINDIR)/private.o $(LINSYS)/scs_matrix.o $(LINSYS)/csparse.o $(LINSYS)/gpu/gpu.o
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -shared -Wl,$(SONAME),$(@:$(OUT)/%=%) -o $@ $^ $(LDFLAGS) $(CULDFLAGS)

$(OUT)/libscsgpuindir.a: $(SCS_INDIR_O) $(SCS_OBJECTS) $(GPUINDIR)/private.o $(LINSYS)/scs_matrix.o $(LINSYS)/csparse.o $(LINSYS)/gpu/gpu.o
	mkdir -p $(OUT)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

# $(OUT)/demo_socp_gpu_direct: test/random_socp_prob.c $(OUT)/libscsgpudir.a
# 	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CULDFLAGS)

$(OUT)/demo_socp_gpu_indirect: test/random_socp_prob.c $(OUT)/libscsgpuindir.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CULDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(SCS_O) $(SCS_INDIR_O) $(SCS_OBJECTS) $(AMD_OBJS) $(LDL_OBJS) $(LINSYS)/*.o $(DIRSRC)/*.o $(INDIRSRC)/*.o $(GPUDIR)/*.o $(GPUINDIR)/*.o $(LINSYS)/gpu/*.o
	@rm -rf $(OUT)/*.dSYM
	@rm -rf matlab/*.mex*
	@rm -rf .idea
	@rm -rf python/*.pyc
	@rm -rf python/build
purge: clean
	@rm -rf $(OUT)

INSTALL_INC_FILES = $(INC_FILES)

INSTALL_TARGETS = $(OUT)/libscsdir.a $(OUT)/libscsindir.a $(OUT)/libscsdir.$(SHARED) $(OUT)/libscsindir.$(SHARED)
INSTALL_GPU_TARGETS = $(OUT)/libscsgpuindir.a $(OUT)/libscsgpuindir.$(SHARED) # $(OUT)/libscsgpudir.a $(OUT)/libscsgpudir.$(SHARED)

INSTALL_INC_DIR = $(DESTDIR)$(PREFIX)/include/scs/
INSTALL_LIB_DIR = $(DESTDIR)$(PREFIX)/lib/

.PHONY: install install_gpu direct indirect
install: $(INSTALL_INC_FILES) $(INSTALL_TARGETS)
	$(INSTALL) -d $(INSTALL_INC_DIR) $(INSTALL_LIB_DIR)
	$(INSTALL) -m 644 $(INSTALL_INC_FILES) $(INSTALL_INC_DIR)
	$(INSTALL) -m 644 $(INSTALL_TARGETS) $(INSTALL_LIB_DIR)
install_gpu: $(INSTALL_INC_FILES) $(INSTALL_GPU_TARGETS)
	$(INSTALL) -d $(INSTALL_INC_DIR) $(INSTALL_LIB_DIR)
	$(INSTALL) -m 644 $(INSTALL_INC_FILES) $(INSTALL_INC_DIR)
	$(INSTALL) -m 644 $(INSTALL_GPU_TARGETS) $(INSTALL_LIB_DIR)
direct:$(OUT)/libscsdir.$(SHARED) $(OUT)/demo_socp_direct $(OUT)/run_from_file_direct $(OUT)/run_tests_direct
indirect:$(OUT)/libscsindir.$(SHARED) $(OUT)/demo_socp_indirect $(OUT)/run_from_file_indirect $(OUT)/run_tests_indirect
