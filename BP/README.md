original_BP will randomly generate some error vectors and decode them by belief propagation.    

The executable programs cla_myparallel, cla_real_parallel, cla_sequential use different BP schedules. To run these program, use    

c
 
An example for the parity check matrix file:
  
3 2  
1 2   
2 3  
Which is the matrix:  
  1 1 0  
  0 1 1  
  
The first line is n n-k, then is the sparse matrix form. The first element is labeled 1 rather than 0.
M2, M3... are some matrices from MayKay's website.
  
  
BP_analysis.cpp will output all the messages after each iteration.  
  bptest is the executable program.  
  `./bptest  <file_name which stored a parity check matrix>  <error probability>  <number of codewords>  <max  iterations>`

  


