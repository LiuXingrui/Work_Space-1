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
range_list=[0.9]
r=0.9
p_range=[0.1,0.01,0.001]
d_range=range(3,13)
for p in p_range:
    data_file_list_s=[]
    data_file_list_p=[]
    for d in d_range:
        n=d*d*2
        data_file_list_s.append("n%d_p%.3f_serial_2.data"%(n,p))
        data_file_list_p.append("n%d_p%.3f_parallel_2.data"%(n,p))
        
    combine( data_file_list_s,"n_p%.3f_serial_2.data"%(p))
    combine( data_file_list_p,"n_p%.3f_parallel_2.data"%(p))
