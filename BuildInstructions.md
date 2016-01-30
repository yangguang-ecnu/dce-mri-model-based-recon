# Introduction #

To use this code, you need to have CUDA and MATLAB installed and running. We are using and testing against CUDA 3.1 and MATLAB R2010a on linux machines, but other versions and operating systems may work as well.

We are also not using the standard CUDA CMake build system, so at this point in time you'll have to manually edit the makefile to point the include and lib flags to the appropriate places for your system.


# Details #

Steps to run:
  1. Download code.
  1. Update makefile flags CUDA\_INCLUDE and CUDA\_LIBRARIES to be accurate for your system.
  1. Run make ( should produce a .so and a .mexa64 file )
  1. Copy the .so to somewhere on your library path and/or update LD\_LIBRARY\_PATH to point at the current dir.
  1. From MATLAB, you should now have access to the dce\_mri\_kernel( ... ) function.

This functions takes a number of arguments and returns a 3D matrix of complex floats. Details of the input and output are outlined in the header comments of dce\_mri\_mex.c.