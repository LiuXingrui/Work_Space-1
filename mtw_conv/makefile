WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 


conv: read_mtx.cpp
	${CPP}  `pkg-config --cflags itpp` -o conv.out read_mtx.cpp  `pkg-config --libs itpp`
