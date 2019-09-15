all: z_density.cpp triangular_kernel_polynomial.o
	g++ -g z_density.cpp triangular_kernel_polynomial.o -lpng -lGL -lGLU -lglut -o z_density

# Get Multivariate Normally Distributed Samples
triangular_kernel_polynomial.o: triangular_kernel_polynomial.cpp 
	g++ -g -c triangular_kernel_polynomial.cpp

