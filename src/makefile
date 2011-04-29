#----------------------------------------------------------------------------
# Authors: Ben Felsted, Simon Williams, Cheng Ye
#    Date: Spring 2011
#----------------------------------------------------------------------------

# Reference calls:
# nvcc --shared -o libsimple.so simple.cu --compiler-options '-fPIC' -I /usr/local/apps/cuda/3.1/cuda/include -L /usr/local/apps/cuda/3.1/cuda/lib64 -lcudart
# mex mlsimple2.c -L. -lsimple

CC=/usr/bin/gcc

# setup CUDA flags for shared lib generation
CUDA_CC        = nvcc
CUDA_FLAGS     = --shared --compiler-options '-v -fPIC -O3' --compiler-bindir /usr/bin
CUDA_INCLUDE   = -I/usr/local/apps/cuda/3.1/cuda/include -I/usr/local/apps/cuda/SDK2/C/common/inc
CUDA_LIBRARIES = -L/usr/local/apps/cuda/3.1/cuda/lib64 -L/usr/local/apps/cuda/SDK2/C/lib
CUDA_LIBS      = -lcudart -lcutil
CUDA_SRC		   = dce_mri_kernel.cu

# setup mex flags for mex function generation
MEX_CC  = mex
MEX_SRC = dce_mri_mex.c

all: kernel mex

kernel:
	$(CUDA_CC) -o libdcemri.so $(CUDA_SRC) $(CUDA_FLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARIES) $(CUDA_LIBS)

mex: kernel
	$(MEX_CC) $(MEX_SRC) -L. -ldcemri

clean:
	rm -rf *.o *.so *.mexa64
