# MAKEFILE for scs java libs
ROOT = ..
include $(ROOT)/scs.mk

JAVA_SRC = src
BIN = bin
OBJECTS = $(ROOT)/src/scs.o $(ROOT)/src/util.o $(ROOT)/src/cones.o $(ROOT)/src/cs.o $(ROOT)/src/linAlg.o $(ROOT)/src/ctrlc.o $(ROOT)/src/scs_version.o $(ROOT)/$(LINSYS)/common.o

AMD_SOURCE = $(wildcard $(ROOT)/$(DIRSRCEXT)/amd_*.c)
DIRECT_OBJECTS = $(ROOT)/$(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o) $(ROOT)/$(DIRSRC)/private.o
INDIRECT_OBJECTS = $(ROOT)/$(INDIRSRC)/private.o

.PHONY: default

default: scs $(BIN)/libjscsdir.$(SHARED) $(BIN)/libjscsindir.$(SHARED) scs.jar
	@echo "**********************************************************************************"
	@echo "Successfully compiled scs java, copyright Brendan O'Donoghue 2012."
	@echo "To test, type 'make testproblem'."
	@echo "To use in your own java projects add scs.jar to the classpath, and set"
	@echo " -Djava.library.path=<your location libraries libjscsdir and libjscsindir>"
	@echo "**********************************************************************************"

scs_purge:
	( cd $(ROOT) && $(MAKE) purge )

scs:
	( cd $(ROOT) && $(MAKE) )

$(BIN)/libjscsdir.$(SHARED): $(OBJECTS) $(DIRECT_OBJECTS) scsjnidir.o
	mkdir -p $(BIN)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

scsjnidir.o:
	$(CC) -c -I$(ROOT)/include -I$(ROOT) -fPIC scsjni.c -o $@

$(BIN)/libjscsindir.$(SHARED): $(OBJECTS) $(INDIRECT_OBJECTS) scsjniindir.o
	mkdir -p $(BIN)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

scsjniindir.o:
	$(CC) -c -DINDIRECTJ -I$(ROOT)/include -I$(ROOT) -fPIC scsjni.c -o $@

scs.jar:
	javac $(JAVA_SRC)/org/scs/*.java
	jar -cf scs.jar -C $(JAVA_SRC) org

testproblem: default
	javac -cp scs.jar $(JAVA_SRC)/test/*.java
	java -cp scs.jar:$(JAVA_SRC)/test -Djava.library.path=$(BIN) ScsTestProblem

.PHONY: clean purge
clean:
	rm -rf scs.jar $(BIN) *.o $(JAVA_SRC)/org/scs/*.class $(JAVA_SRC)/test/*.class
purge: scs_purge clean

