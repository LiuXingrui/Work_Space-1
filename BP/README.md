BP_basic_func1.cpp and BP_quantum_func1.cpp contains the functions used in the programs. original_BP.cpp is for classical BP, quantum_BP_x_only.cpp is for quantum BP. We can ignore most other files. For compiling, check Makefile, for compiling in cluster, use   
hpcc1= ``` `itpp-config --cflag\s` ``` 
hpcc2= ``` `itpp-config --libs ` ```    
instead of  
``` `pkg-config --cflags itpp` ``` and   
``` `pkg-config --libs itpp` ```   

original_BP will randomly generate some error vectors and decode them by belief propagation.    

The executable programs cla_myparallel, cla_real_parallel, cla_sequential....(Some of them are old editions, it is better recompile before running) use different BP schedules. To run these program, use    

  `./new_sequential  <H_file> <data_file> <pmin> <pmax>   <number of codewords>  <max  iterations>`   
  The data stored in data_file are 0:n 1:fail_rate 2:avg_p  3:avg_iter 4:number_of_suc_decoding 5:bit_error_rate_after_decoding
  
 
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
  
  
  xqBP will decode CSS code with error probability distributed between pmin and pmax:  
  `./xqBP <Hx_file> <Hz_file>  <pavg/wt> <range> <max number of failed decoding/num_of_cws><lmax> <data_file><debug><channel><alpha><decode_p><decode_prange><lambda>`    
  if <max number of failed decoding/num_of_cws> is less than 0, then that is num_of_cws. If it is greater than 0, that is this max_num_of_failed_decodings, when the program produces such number of failed decodings, it stopped.  
  if <pavg/wt> is greater or equal than 1, it is explained as the weight of the errors and pavg will be setted to be wt/n, for example, if  <pavg/wt> =2, the program will only generate wt=2 errors.   
  <channel> :
   channel==0: bsc, only decode z-errors(use Hx)  
   channel==1: depolarizing, pavg is the depolarizing rate  
  it is better not to use channel==1, xqBP cannot deal with it.  
  
  <debug>:  
  debug%2==1 : if reach maximum iteratios, try another p_dis, decode it again.  
  (debug/2)%2==1 : parallel  schedule  
   (debug/4)%2==1: print the real_e, output_e when reach max iterations  
   (debug/8)%2==1: use OSD after BP fails.  
   (debug/16)%2==1: use LLR_avg.   (don't use this now, it won't work.)
   0001=1: try again  
   0010=2: para  
   0011=3: try again and para  
   0100=4: print mes  
   0101=5: print and try again  
   0110=6: print and  para
   0111=7: print and para and try again  
   1001=9: OSD and try again  
   1011=11: para and OSD and try again  

   <alpha> is for testing another method.  
   <lambda> won't work now, just set it=1 is ok.
     
   decode_p/decode_prange are the p/range for decoding, pavg and range are for the real p. In fact, it seems the real p_range does not matter???
   
   
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
   12: decode_p  
   13:decode_prange    
   14:alpha  
   15: OSD_suc  
   16:lambda  
   

  gen_cyclic.cpp will create a rank=r-1 parity check matrix file for a cyclic code, r is the number of rows:  
  `./gen_cyc <file for storing> <n> <h_k> <h_k-1> ...<h_0>`


 HT.cpp will write the transposed H to a file:  
 `./tran <H> <H^T>`
 
 `./display <H>` prints a dense mat.  

 Check create_toric_mat.py and Makefile for more commands.
 
 pyplot is for ploting.  
 
