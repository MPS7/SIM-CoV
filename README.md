SIM-CoV v1.1b
=============

SIM-CoV is a hybrid multi-scale simulator for COVID-19 transmission dynamics.


# Dependencies

The SIM-CoV code was developed and tested under Ubuntu 16.04. There exist few dependencies that need to be installed on the OS before compiling and running SIM-CoV-2020:

- GCC: http://gcc.gnu.org/
	- apt-get install gcc

- GNU Scientific Library (GSL) https://www.gnu.org/software/gsl/
	- download and install the GSL package: ftp://ftp.gnu.org/gnu/gsl/

# Compiling and running

Use the following commands to compile, build and execute:

g++ -I/usr/local/include -c main.cc 
g++ -L/usr/local/lib main.o -lgsl -lgslcblas -lm

./a.out


# Post-processing and visualization

VTK files are generated in the vtk_particle and vtk_PDE folders upon running the program. We recommand using ParaView to properly visualize the simulations. The generated file data.txt contains the number of the different agent types over time. The file data2.txt contains demographic characteristics of the sampled population.
