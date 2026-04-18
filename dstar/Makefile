CXXFLAGS  += -I../include -I../include/star \
	     -D_SRC_='"no_source_available"' \
	     -DHAVE_CONFIG_H  -DTOOLBOX 
LDLIBS	  += \
	     -L. -ldstar \
	     -L../sstar -lsstar \
	     -L../node -lnode \
	     -L../node/dyn -ldyn \
	     -L../std -lstd  \
	     -lm

CFLAGS    += -O
# OpenMP: enabled on Linux with GCC (-fopenmp).
# Apple Clang on macOS does not support -fopenmp without libomp;
# the binary falls back to single-threaded automatically via #ifdef _OPENMP guards.
# On the cluster (Linux/GCC) this gives near-linear speedup with core count.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS  += -O -fopenmp
    LDLIBS    += -fopenmp
else
    CXXFLAGS  += -O
endif

EXE	= SeBa Normalisation
DIRS	= init io starclass stardyn util

all:	lib $(EXE)

Makefile.inc:	Makefile.inc.conf
	@sed s#__BASE__#$(PWD)#g Makefile.inc.conf >Makefile.inc

lib:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done

clean:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done
	/bin/rm -f *~ *.a $(EXE) Makefile.inc

