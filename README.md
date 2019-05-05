# SincDVR

This program is aimed to run simulations of a quantum system in one dimension. It is based on a sinc function basis.

Currently this program is able to read arbitrary combination of elementary functions (by multiplying or adding them) as the initial wave function or potential function. The spaces of the grids, the number of the grid points as well as the time interval needs to be modified manually. The output only records the values of the wave function on the grid points at each step, which means the mean value of any observable needs to be calculated manually.

The program 'seems' to be written in C++, but actually it's completely C style (yes it is). 

The program greatly utilizes GSL (GNU Scientific Library, https://www.gnu.org/software/gsl/). Thus you need to first install the whole package of this library to compile & run the program. Currently the compiling process is written as a shell script stored in /build.
