import os
d_list=[5,7,9]
for i in range(5,11,1):
    
    alpha=0.1*i
    command= "rm alpha%.1f_serial.data "%(alpha)
    os.system(command)
    command= "rm alpha%.1f_parallel.data "%(alpha)
    os.system(command)
    
    for d in d_list:
        command= "rm alpha%.1f_d%d_serial.data "%(alpha,d)
        os.system(command)
        command= "rm alpha%.1f_d%d_parallel.data "%(alpha,d)
        os.system(command)
        
for i in range(5,11,1):
    alpha=0.1*i
    for d in d_list:    
        command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.05 0.9 10000 20  alpha%.1f_serial.data 1 0 %f 0.05 0.9 "%(d,d,alpha,alpha)
        os.system(command) 
        command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.05 0.9 10000 20  alpha%.1f_parallel.data 1 0 %f 0.05 0.9 "%(d,d,alpha,alpha)
        os.system(command)  
