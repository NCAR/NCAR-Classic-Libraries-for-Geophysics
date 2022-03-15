
include make.inc

all: libspherepack testspherepack

libspherepack:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) )

testspherepack:
	( cd ./test; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )
