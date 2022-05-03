WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 


dist: dist.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp distance.cpp lib.cpp mm_read.cpp mm_write.cpp mmio.c
	${CPP}  `pkg-config --cflags itpp` -o distance dist.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp distance.cpp lib.cpp mm_read.cpp mm_write.cpp mmio.c `pkg-config --libs itpp`
