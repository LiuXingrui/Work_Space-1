
d_list=range(5,7,9)

for alpha in range(0.5,1.1,0.1):
    for d in d_list:    
        command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.01 0.9 10000 100  alpha%.1f_d%d_serial.data 1 0 %f 0.01 0.9 "%(d,d,alpha,d,alpha)
        os.system(command) 
        command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.01 0.9 10000 100  alpha%.1f_d%d_parallel.data 1 0 %f 0.01 0.9 "%(d,d,alpha,d,alpha)
        os.system(command)  
