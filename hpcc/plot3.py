import os
index_list=range(1,6)
n_list=[10,58,106,166,202]

pavg_list = [0.1,0.01,0.001]
range_list=[0,0.9]

alpha=1

for i in range(1,6):
    n=n_list[i-1]
    
    for p in pavg_list:
        for r in range_list:
            
            command="./xqBP GB_w4_%d_X   GB_w4_%d_Z  %f %f 20000 100  bicycle_n%d_p%.3f_range%.1f.data 1 0 1 %f %f"%(i,i,p,r,n,p,r,p,r)
            os.system(command) 
      
