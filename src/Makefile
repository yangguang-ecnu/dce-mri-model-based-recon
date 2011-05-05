#----------------------------------------------------------------------------
# Authors: Ben Felsted, Simon Williams, Cheng Ye
#    Date: Spring 2011
#----------------------------------------------------------------------------

# Default targets to make
all: kernel mex


# If you want to change make variables from their defaults, edit 
# 'Makefile.in' or pass arguments as assignments to the make program. The 
# latter approach can be done as the following example shows:
.PHONY : example 
example:  
	make CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda/3.1/cuda \
		 CUDA_SDK_ROOT_DIR=/usr/local/cuda/SDK2/C \
		 VERBOSE=-v
# This example won't work in the CADE lab btw

#-----------------------------------------------------------------------------
# Source all of the make variable parameters into this file
include Makefile.in

# Reference calls:
# $ nvcc --shared -o libsimple.so simple.cu --compiler-options '-fPIC' -I /usr/local/apps/cuda/3.1/cuda/include -L /usr/local/apps/cuda/3.1/cuda/lib64 -lcudart
# $ mex mlsimple2.c -L. -lsimple -Wl,-rpath .

kernel: 
	$(CUDA_CC) -o libdcemri.so $(CUDA_SRC) $(CUDA_FLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARY) $(CUDA_LIBS)
	
mex: kernel
	$(MEX_CC) $(MEX_SRC) $(MEX_LIBRARY) $(MEX_LIBS) $(MEX_CC_FLAGS)

clean:
	rm -rf *.o *.so *.mexa64


# ------------ A few relevent gcc compiler options ------------------# {{{1
#            -ffast-math                                             # {{{2
#                Sets -fno-math-errno, -funsafe-math-optimizations, -ffi-
#                nite-math-only, -fno-rounding-math, -fno-signaling-nans
#                and -fcx-limited-range.
#
#                This option causes the preprocessor macro
#                "__FAST_MATH__" to be defined.
#
#                This option is not turned on by any -O option since it
#                can result in incorrect output for programs which depend
#                on an exact implementation of IEEE or ISO rules/specifi-
#                cations for math functions. It may, however, yield
#                faster code for programs that do not require the guaran-
#                tees of these specifications.
#
#
#               i ---------------------------------------------------# }}}2
#            -funsafe-math-optimizations                             # {{{2
#                Allow optimizations for floating-point arithmetic that
#                (a) assume that arguments and results are valid and (b)
#                may violate IEEE or ANSI standards.  When used at
#                link-time, it may include libraries or startup files
#                that change the default FPU control word or other simi-
#                lar optimizations.
#
#                This option is not turned on by any -O option since it
#                can result in incorrect output for programs which depend
#                on an exact implementation of IEEE or ISO rules/specifi-
#                cations for math functions. It may, however, yield
#                faster code for programs that do not require the guaran-
#                tees of these specifications.  Enables
#                -fno-signed-zeros, -fno-trapping-math, -fassocia-
#                tive-math and -freciprocal-math.
#
#                The default is -fno-unsafe-math-optimizations.
#
#
#                ----------------------------------------------------# }}}2
#            -fPIC                                                   # {{{2
#                If supported for the target machine, emit position-
#                independent code, suitable for dynamic linking and
#                avoiding any limit on the size of the global offset
#                table.  This option makes a difference on the m68k,
#                PowerPC and SPARC.
#
#                Position-independent code requires special support,
#                and therefore works only on certain machines.
#
#
#                ---------------------------------------------------# }}}2
# }}}1
# vim:fdm=marker
