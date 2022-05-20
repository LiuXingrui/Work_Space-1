import os
wt_list=[3,5,6,10]
n_list=[18,50,72,200]
range_list=[0,0.1,0.2,0.3,0.6,0.9,1.0]

for max_wt in wt_list:
    n=max_wt*max_wt*2  
    wt_range=range(1,max_wt+1)
    for r in range_list:
        os.system("rm wt_n%d_%f.txt"%(n,r))
        os.system("rm wt_n%d_%.1f.txt"%(n,r))

for max_wt in wt_list:
    n=max_wt*max_wt*2  
    wt_range=range(1,max_wt+1)
    
    for wt in wt_range:
        for r in range_list:
            command="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 %d %f 1000 20 wt_n%d_%.1f.txt "%(max_wt,max_wt,wt,r,n,r)
            os.system(command)  
      
    
