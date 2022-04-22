test: test.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	g++   `pkg-config --cflags itpp` -o test test.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`
 

gen_Hx:gen_surface.cpp
	g++   -o gen_surface_star gen_surface.cpp

hypergraph: gen_HPC.cpp BP_basic_func1.cpp
	g++   `pkg-config --cflags itpp` -o HPC gen_HPC.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`
 
qtest1: quantum_test1.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	g++   `pkg-config --cflags itpp` -o qtest1 quantum_test1.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

cla_BP:BP_basic_func1.cpp original_BP.cpp BP_quantum_func1.cpp
	g++   `pkg-config --cflags itpp` -o new_sequential original_BP.cpp  BP_basic_func1.cpp BP_quantum_func1.cpp  `pkg-config --libs itpp`

quan_BP:quantum_BP.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
		g++   `pkg-config --cflags itpp` -o qBP quantum_BP.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

gen_cyc: gen_cyclic.cpp  BP_basic_func1.cpp
	g++   `pkg-config --cflags itpp` -o gen_cyc gen_cyclic.cpp  BP_basic_func1.cpp  `pkg-config --libs itpp`

tran: HT.cpp  BP_basic_func1.cpp
	g++   `pkg-config --cflags itpp` -o tran HT.cpp  BP_basic_func1.cpp  `pkg-config --libs itpp`
