import numpy as np
import pandas as pd
def combine(data_file_list,target_file):
    
    n_data=len(data_file_list)
    newdata=np.zeros((n_data,14))
    
    same=[0,1,3,4,12,13]
    add=[6,7,8,9]
    divide=[2,10,11] #2=(7-6)/7, 10=8/7, 11=9/7
    avg=5
    row_index=0
    
    for data_file in data_file_list:

        data=pd.read_csv(data_file, sep='\s+',header=None)
        data=data.to_numpy()
        rows=data.shape[0]    

        for i in same:
            newdata[row_index,i]=data[0,i]

        for i in range(rows):
            newdata[row_index,5]=newdata[row_index,5]+data[i,5]*data[i,7]
        for j in add:
            newdata[row_index,j]=newdata[row_index,j]+data[i,j]
        
        newdata[row_index,5]=newdata[row_index,5]/float(newdata[row_index,7])
        newdata[row_index,2]=(newdata[row_index,7]-newdata[row_index,6])/float(newdata[row_index,7])
        newdata[row_index,10]=newdata[row_index,8]/float(newdata[row_index,7])
        newdata[row_index,11]=newdata[row_index,9]/float(newdata[row_index,7])
        row_index=row_index+1
    print(newdata)
    np.savetxt(target_file, newdata, fmt='%.6f',delimiter=' ', newline='\n', encoding=None)

#pavg_list = [float(x)*0.001 for x in range(1,11)]
range_list=[0,0.9]
wt_list=[4,6]
p_range=[0.1,0.01,0.001]
n_range=[10,58,106,166,202]
for p in p_range:
    for wt in wt_list:
        for r in range_list:
            
            data_file_list=[]
 
            for i in range(1,6):
                n=n_range[i-1]
                data_file_list.append("bicycle_w%d_n%d_p%.3f_range%.1f.data"%(wt,n,p,r))
       
            combine( data_file_list,"bicycle_w%d_p%.3f_range%.1f.data"%(wt,p,r))
  
