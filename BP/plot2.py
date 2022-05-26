import os
wt_list=[17]
n_list=[578]
pavg_list = [0.025,0.05]
range_list=[0.9]
alpha=0.65

for max_wt in wt_list:
    n=max_wt*max_wt*2  
    for r in range_list:
        os.system("rm 578_%.1f_s2_alpha.data"%r)
        os.system("rm 578_%.1f_s2.data"%r)
        os.system("rm 578_%.1f_p2.data"%r)


for max_wt in wt_list:
    n=max_wt*max_wt*2  
    
    for p in pavg_list:
        for r in range_list:
            command="./xqBP toric_Hx_n17_cyc11 toric_Hz_n17_cyc11 %f %f 1000 100  578_%.1f_s2_alpha.data 1 0 %f "%(p,r,r,alpha)
            os.system(command) 
            command="./xqBP toric_Hx_n17_cyc11 toric_Hz_n17_cyc11 %f %f 1000 100  578_%.1f_s2.data 1 0 1 "%(p,r,r)
            os.system(command) 
   
      
    
