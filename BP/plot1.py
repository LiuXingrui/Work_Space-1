import os
wt_list=[10]
n_list=[200]
pavg_list = [float(x)*0.001 for x in range(1,11)]
range_list=[0,0.1,0.2,0.3,0.6,0.9] 


for max_wt in wt_list:
    n=max_wt*max_wt*2  
    
    for p in pavg_list:
        for r in range_list:
            command="./xqBP toric_Hx_n10_cyc11 toric_Hz_n10_cyc11 %f %f 20 100  n200_range%.1f_p%.3f_serial2.data 1 0 1 %f %f "%(p,r,r,p,p,r)
            os.system(command) 
            command="./xqBP toric_Hx_n10_cyc11 toric_Hz_n10_cyc11 %f %f 20 100  n200_range%.1f_p%.3f_parallel2.data 3 0 1 %f %f"%(p,r,r,p,p,r)
            os.system(command)  
      
    

