CUDA_LIB :=/usr/local/cuda/lib64 -lcuda -lcudart

All: A B C
	g++ -fopenmp FDMLR91.o memory.o forward.o -o test.run -L $(CUDA_LIB)
A:
	nvcc memory.cu  --ptxas-options=-v -c
B:
	nvcc forward.cu  --ptxas-options=-v -c
C:
	g++ FDMLR91.cpp -c -fopenmp
clean:
	rm *.o *.txt *.run
