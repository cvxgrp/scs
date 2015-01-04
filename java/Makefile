# MAKEFILE for scs java libs
ROOT = ..
include $(ROOT)/scs.mk

JAVA_OUT = src
OBJECTS = $(ROOT)/src/scs.o $(ROOT)/src/util.o $(ROOT)/src/cones.o $(ROOT)/src/cs.o $(ROOT)/src/linAlg.o $(ROOT)/src/ctrlc.o $(ROOT)/$(LINSYS)/common.o

AMD_SOURCE = $(wildcard $(ROOT)/$(DIRSRCEXT)/amd_*.c)
DIRECT_OBJECTS = $(ROOT)/$(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o) $(ROOT)/$(DIRSRC)/private.o
INDIRECT_OBJECTS = $(ROOT)/$(INDIRSRC)/private.o

.PHONY: default

default: scs $(JAVA_OUT)/libjscsdir.$(SHARED) $(JAVA_OUT)/libjscsindir.$(SHARED)

scs_purge:
	( cd $(ROOT) && $(MAKE) purge )

scs:
	( cd $(ROOT) && $(MAKE) )

$(JAVA_OUT)/libjscsdir.$(SHARED): $(OBJECTS) $(DIRECT_OBJECTS) scsjnidir.o
	$(CC) -shared -o $@ $^ $(LDFLAGS)

scsjnidir.o:
	$(CC) -c -I$(ROOT)/include -I$(ROOT) scsjni.c -o $@

$(JAVA_OUT)/libjscsindir.$(SHARED): $(OBJECTS) $(INDIRECT_OBJECTS) scsjniindir.o
	$(CC) -shared -o $@ $^ $(LDFLAGS)

scsjniindir.o:
	$(CC) -c -DINDIRECTJ -I$(ROOT)/include -I$(ROOT) scsjni.c -o $@

testproblem: default
	javac -cp $(JAVA_OUT) $(JAVA_OUT)/scs/TestProblem.java
	java -cp $(JAVA_OUT) -Djava.library.path=$(JAVA_OUT) scs.TestProblem

.PHONY: clean purge
clean:
	rm -rf $(JAVA_OUT)/libjscsdir.$(SHARED) $(JAVA_OUT)/libjscsindir.$(SHARED) *.o $(JAVA_OUT)/scs/*.class
purge: scs_purge clean

