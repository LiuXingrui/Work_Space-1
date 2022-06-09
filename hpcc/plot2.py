import os
d_list=range(3,13)

pavg_list = [0.1,0.01,0.001]
range_list=[0.9]
r=0.9
alpha=1

for d in d_list:
    n=d*d*2
    
    for p in pavg_list:
            command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 %f 0.9 20000 100  n%d_p%.3f_serial_2.data 1 0 1 %f 0.9"%(d,d,p,n,p,p)
            os.system(command) 
            command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 %f 0.9 20000 100  n%d_p%.3f_parallel_2.data 3 0 1 %f 0.9"%(d,d,p,n,p,p)
            os.system(command) 
   
      
    
