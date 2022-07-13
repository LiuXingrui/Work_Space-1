import os
wt_list=[10]
n_list=[200]
pavg_list = [0.01,0.03,0.07,0.1]
lambda_list=[0.5,0.7,0.9,1] 


for max_wt in wt_list:
    n=max_wt*max_wt*2  
    
    for p in pavg_list:
        for lamb in lambda_list:
            command="./xqBP toric_Hx_n10_cyc11 toric_Hz_n10_cyc11 %f 0.9 1000 100  ./data/lambdatestp%.2f.data 9 0 1 %f 0.9 %f"%(p,p,p,lamb)
            os.system(command) 
    
    
