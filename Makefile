load_all: load.sh
	source load.sh
g+++: circulant_QHP.cpp
	g++ -std=c++11 -o circulant_QHP.out circulant_QHP.cpp
