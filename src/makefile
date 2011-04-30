#----------------------------------------------------------------------------
# Authors: Ben Felsted, Simon Williams, Cheng Ye
#    Date: Spring 2011
#----------------------------------------------------------------------------

# Reference calls:
# nvcc --shared -o libsimple.so simple.cu --compiler-options '-fPIC' -I /usr/local/apps/cuda/3.1/cuda/include -L /usr/local/apps/cuda/3.1/cuda/lib64 -lcudart
# mex mlsimple2.c -L. -lsimple



# --------------------A few gcc compiler options ------------------# {{{1
#       -ffast-math
#           Sets -fno-math-errno, -funsafe-math-optimizations, -ffi-
#           nite-math-only, -fno-rounding-math, -fno-signaling-nans
#           and -fcx-limited-range.
#
#           This option causes the preprocessor macro
#           "__FAST_MATH__" to be defined.
#
#           This option is not turned on by any -O option since it
#           can result in incorrect output for programs which depend
#           on an exact implementation of IEEE or ISO rules/specifi-
#           cations for math functions. It may, however, yield
#           faster code for programs that do not require the guaran-
#           tees of these specifications.
#
#
#       -funsafe-math-optimizations
#           Allow optimizations for floating-point arithmetic that
#           (a) assume that arguments and results are valid and (b)
#           may violate IEEE or ANSI standards.  When used at
#           link-time, it may include libraries or startup files
#           that change the default FPU control word or other simi-
#           lar optimizations.
#
#           This option is not turned on by any -O option since it
#           can result in incorrect output for programs which depend
#           on an exact implementation of IEEE or ISO rules/specifi-
#           cations for math functions. It may, however, yield
#           faster code for programs that do not require the guaran-
#           tees of these specifications.  Enables
#           -fno-signed-zeros, -fno-trapping-math, -fassocia-
#           tive-math and -freciprocal-math.
#
#           The default is -fno-unsafe-math-optimizations.
# }}}1


#CC=/usr/bin/gcc

#VERBOSE       ?= -v
CC_FLAGS       = -fPIC $(VERBOSE)
CC_BINDIR      = /usr/bin
#CC_BINDIR      = /home/felsted/root/usr/bin
 #--linker-options "-rpath ."

# setup CUDA flags for shared lib generation
CUDA_CC        = nvcc
CUDA_FLAGS     = --shared --compiler-options $(CC_FLAGS) --compiler-bindir $(CC_BINDIR) 
CUDA_INCLUDE   = -I/usr/local/apps/cuda/3.1/cuda/include -I/usr/local/apps/cuda/SDK2/C/common/inc
CUDA_LIBRARIES = -L/usr/local/apps/cuda/3.1/cuda/lib64 -L/usr/local/apps/cuda/SDK2/C/lib
CUDA_LIBS      = -lcudart -lcutil
CUDA_SRC       = dce_mri_kernel.cu

# setup mex flags for mex function generation
MEX_CC  = mex
MEX_SRC = dce_mri_mex.c

all: kernel mex

#dce-mri-kernel.cu:
#	$(CUDA_CC) -o dce-mri-kernel.o $(CUDA_SRC) $(CUDA_FLAGS) $(CUDA_INCLUDE) #$(CUDA_LIBRARIES) $(CUDA_LIBS)
# dce-mri-kernel.cu
#kernel: 
#	time $(CUDA_CC) -o libdcemri.so $(CUDA_SRC) $(CUDA_FLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARIES) $(CUDA_LIBS) | cat

kernel: 
	$(CUDA_CC) -o libdcemri.so $(CUDA_SRC) $(CUDA_FLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARIES) $(CUDA_LIBS)
	
mex: kernel
	$(MEX_CC) $(MEX_SRC) -L. -ldcemri -Wl,-rpath . CC=$(CC_BINDIR)/gcc CXX=$(CC_BINDIR)/g++ #-v -Wl,--trace 

clean:
	rm -rf *.o *.so *.mexa64

