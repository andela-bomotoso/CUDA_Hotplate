p6: hotplate.cu
	/usr/local/cuda/bin/nvcc -o p6 hotplate.cu
	
clean:
	rm -f p6 *.o
