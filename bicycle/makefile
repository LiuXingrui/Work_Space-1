WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 

bike: func.cpp func.hpp main.cpp
	${CPP}  `pkg-config --cflags itpp` -o bike.out func.cpp main.cpp  `pkg-config --libs itpp`
