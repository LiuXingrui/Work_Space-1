To compile, use 
g++   `pkg-config --cflags itpp` -o prog_name original_BP.cpp  BP_basic_func1.cpp  `pkg-config --libs itpp`

cla_myparallel, cla_real_parallel, cla_sequential use different BP schedules. To run these program, use 
./cla_parallel  <file stored a parity check matrix> <error probability> <number of codewords> <max  iterations>
An example for the parity check matrix file:
  
3 2
1 2 
2 3
 Which is the matrix:
  1 1 0
  0 1 1
  
The first line is n n-k, then is the sparse matrix form. The first element is labeled 1 rather than 0.
  
 For more details, see 
  https://www.overleaf.com/read/gkwhwsdwrjsj
  The notes are not finished yet.
  
