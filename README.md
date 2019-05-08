# SincDVR

This program is aimed to run simulations of a quantum system in one dimension. It is based on a sinc function basis.

Currently this program is able to read arbitrary combination of elementary functions (by multiplying or adding them) as the initial wave function or potential function. The spaces of the grids, the number of the grid points as well as the time interval needs to be modified manually. 

The program 'seems' to be written in C++, but actually it's completely C style (yes it is). 

The program greatly utilizes GSL (GNU Scientific Library, https://www.gnu.org/software/gsl/). Thus you need to first install the whole package of this library to compile & run the program. Currently the compiling process is written as a shell script stored in /build.

Two examples of writing an input file are given as txt files in /input. Use flag "-i" along with the path & file name of the input file to make an input, namely

$ ./DVR -i input/input.txt

Default output file is stored as "NONAME.txt" in /output directory. Use flag "-o" to specify the output file name or directory. Tails of ".swp" are automatically added in the output file if the file name for the output already exists, and tail of ".log" for log files containing information of Hamiltonian Matrix and other unnecessary information.
