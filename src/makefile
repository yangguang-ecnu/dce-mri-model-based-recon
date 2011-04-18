#----------------------------------------------------------------------------
# Authors: Ben Felsted, Simon Williams, Cheng Ye
#    Date: Spring 2011
#----------------------------------------------------------------------------

# Reference calls:
# nvcc --shared -o libsimple.so simple.cu --compiler-options '-fPIC' -I /usr/local/apps/cuda/3.1/cuda/include -L /usr/local/apps/cuda/3.1/cuda/lib64 -lcudart
# mex mlsimple2.c -L. -lsimple

# setup CUDA flags for shared lib generation
CUDA_CC        = nvcc
CUDA_FLAGS     = --shared --compiler-options '-fPIC'
CUDA_INCLUDE   = /usr/local/apps/cuda/3.1/cuda/include
CUDA_LIBRARIES = /usr/local/apps/cuda/3.1/cuda/lib64
CUDA_LIBS      = -lcudart
CUDA_SRC		   = dce-mri-kernel.cu
CUDA_OUT       = libdcemri.so

# setup mex flags for mex function generation
MEX_CC  = mex
MEX_SRC = dce-mri-mex.c

all: kernel mex

kernel:
	$(CUDA_CC) -o $(CUDA_OUT) $(CUDA_SRC) $(CUDA_FLAGS) -I $(CUDA_INCLUDE) -L $(CUDA_LIBRARIES) $(CUDA_LIBS)

mex: kernel
	$(MEX_CC) $(MEX_SRC) -L. -l$(CUDA_OUT)

clean:
	rm -rf *.o *.so *.mex64
