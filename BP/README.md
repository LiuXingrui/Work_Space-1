original_BP will randomly generate some error vectors and decode them by belief propagation.    

The executable programs cla_myparallel, cla_real_parallel, cla_sequential use different BP schedules. To run these program, use    

  `./new_sequential  <H_file> <data_file> <pmin> <pmax>   <number of codewords>  <max  iterations>`   
  The data stored in data_file are n fail_rate avg_p  avg_iter number_of_suc_decoding bit_error_rate_after_decoding   
 
An example for the parity check matrix file:
  
3 2  
1 2   
2 3  
Which is the matrix:  
  1 1 0  
  0 1 1  
  
The first line is n n-k, then is the sparse matrix form. The first element is labeled 1 rather than 0.
M2, M3... are some matrices from MayKay's website.
  
  gen_HPC.cpp will generate Hx and Hz for hypergraph product code.    
  `./HPC <file stored H1> <file stored H2> <file stored Hx> <file stored Hz>`
  
  
  qBP will decode CSS code with error probability distributed between pmin and pmax:  
  `./qBP <Hx_file> <Hz_file>  <pavg/wt> <range> <number of codewords><lmax> <data_file><debug><channel>`    
  debug%2==1: if reach maximum iteratios, try another p_dis, decode it again.  
  (debug/2)%2==1: parallel  schedule  
   (debug/4)%2==1: print messages after each iterations if reach max iterations  
   channel==0: bsc, only decode z-errors(use Hx)  
   channel==1: depolarizing, pavg is the depolarizing rate  
   
   
   The data stored in data_file are  n d fail_rate pavg/wt range avg_iter_for_suc num_of_suc_dec num_of_cws syn_fail max_fail syn_fail_rate max_fail_rate:  
   0:n  
   1:d  (have not calculate it in this prog yet, so outputs for d are setted to -1)
   2: fail_rate  
   3:pavg/wt  
   4: range  
   5:avg_iter_for_suc  
   6:num_of_suc_dec  
   7:  num_of_cws  
   8: syn_fail  
   9:max_fail  
   10:syn_fail_rate  
   11: max_fail_rate  
  
   if "pavg/wt">=1, then it will be explained as weight of errors and pavg=wt/n.
   
   qBPx only decode z-errors (use Hx).  
   qBPx2 will decode the same error again use another random p_distribution if reaches maximum iterations.
   
   
   aqBP gives more information.   aqBP_try_again works same way as qBPx2, aqBP_print_mes print messages after each iteration if the decodeing reaches the maximum iterations.

  gen_cyclic.cpp will create a rank=r-1 parity check matrix file for a cyclic code, r is the number of rows:  
  `./gen_cyc <file for storing> <n> <h_k> <h_k-1> ...<h_0>`


 HT.cpp will write the transposed H to a file:  
 `./tran <H> <H^T>`
 
 `./display <H>` prints a dense mat.  

 Check create_toric_mat.py and Makefile for more commands.
 
 pyplot is for ploting.  
 
