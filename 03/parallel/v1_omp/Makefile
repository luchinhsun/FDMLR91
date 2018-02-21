CUDA_LIB :=/usr/local/cuda/lib64 -lcuda -lcudart

comp := g++
omp := -fopenmp
#comp := icc
#omp := -qopenmp
#opmn := -O3

All: A B C
	$(comp) $(opmn) $(omp) FDMLR91.o memory.o forward.o -o test.run -L $(CUDA_LIB)
A:
	nvcc $(opmn) memory.cu --ptxas-options=-v -c
B:
	nvcc $(opmn) forward.cu --ptxas-options=-v -c
C:
	$(comp) $(opmn) $(omp) FDMLR91.cpp -c
clean:
	rm *.o *.txt *.run
