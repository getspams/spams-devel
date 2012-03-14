export ARCH=intel64
export CXX=icpc
export MATH_LIB=mkl
export DEBUG=off
export MAX_THREADS=64
export CXXARCHFLAG=-axSSE3,SSE4.1,SSE4.2,AVX

############################################
# You should not change the following part #
############################################

export DECLARE=-DMAX_THREADS=$(MAX_THREADS) -D_REENTRANT 
export INCLUDE=-I$(MATLABPATH)/extern/include/
export INCLUDE+=-I$(MATLABPATH)/extern/src/

export DECLARE+=-DEM64T
export MEXEXT=mexa64
export MATLABVERSION=glnxa64

export INCLUDE+=-I$(ICPCPATH)/include/ #-I$(GXX_INCLUDE) -I$(GXX_INCLUDE2)
export INCLUDELIB=-L$(ICPCPATH)/lib/intel64/

#export LIBS=-Bstatic -lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-Bdynamic -lmkl_lapack
#export LIBS+=-Bdynamic -liomp5 -lpthread
PATHLIB=$(MKLPATH)/lib/intel64/
export LIBS+= -Wl,--start-group ${PATHLIB}/libmkl_intel_lp64.a ${PATHLIB}/libmkl_sequential.a ${PATHLIB}/libmkl_core.a -Wl,--end-group -liomp5 -lpthread

export INCLUDE+=-I$(MKLPATH)/include/ 
export INCLUDELIB+=-L$(MKLPATH)/lib/$(ARCH)/ -L$(MKLPATH)/../lib/$(ARCH)/ -L./
#export DECLARE+=-DHAVE_MKL

#export LIBS+=-Bstatic -lsvml
export CXXARCHFLAG+= -fPIC -pipe -ipo -w -w0 -O3 -fomit-frame-pointer -no-prec-div -DNDEBUG -openmp -parallel -fno-alias -fno-fnalias -align -falign-functions -fp-model fast -funroll-loops -cxxlib-gcc #-gcc-version=410
export CXX=$(ICPCPATH)/bin/intel64/icpc
export MEXLINK=-shared -Wl,--version-script,$(MATLABPATH)/extern/lib/$(MATLABVERSION)/mexFunction.map
export mexversionsrc=$(MATLABPATH)/extern/src/mexversion.c
export CFLAGS=$(CXXARCHFLAG) $(DECLARE) 
export NAME=-o
export LIBNAME=libSPAMS.a
export VER=linuxmac

# Build rules

headers=../linalg/utils.h ../linalg/mkl_template.h ../linalg/misc.h ../linalg/linalg.h
objects=mexversion.o

all: 
	mkdir -p build
	$(MAKE) toolbox_decomp
	$(MAKE) toolbox_prox
	$(MAKE) toolbox_linalg
	$(MAKE) toolbox_dictLearn

toolbox_linalg:
	@(cd linalg && $(MAKE))

toolbox_prox:
	@(cd prox && $(MAKE))

toolbox_decomp:
	@(cd decomp && $(MAKE))

toolbox_dictLearn:
	@(cd dictLearn && $(MAKE))

toolbox_test:
	@(cd test && $(MAKE))

clean:
	@(cd linalg && $(MAKE) clean)
	@(cd prox && $(MAKE) clean)
	@(cd decomp && $(MAKE) clean)
	@(cd dictLearn && $(MAKE) clean)
	rm build/*.o build/*.mexa64
	rm build/test


