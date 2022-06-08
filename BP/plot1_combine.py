import numpy as np
import pandas as pd
def combine(data_file,target_file):

    data=pd.read_csv(data_file, sep='\s+',header=None)
    data=data.to_numpy()

    rows=data.shape[0]
    newdata=np.zeros((1,14))

    same=[0,1,3,4,12,13]
    add=[6,7,8,9]
    divide=[2,10,11] #2=(7-6)/7, 10=8/7, 11=9/7
    avg=5

    for i in same:
        newdata[0,i]=data[0,i]

    for i in range(rows):
        newdata[0,5]=newdata[0,5]+data[i,5]*data[i,7]
        for j in add:
            newdata[0,j]=newdata[0,j]+data[i,j]
        

    newdata[0,5]=newdata[0,5]/float(newdata[0,7])
    newdata[0,2]=(newdata[0,7]-newdata[0,6])/float(newdata[0,7])
    newdata[0,10]=newdata[0,8]/float(newdata[0,7])
    newdata[0,11]=newdata[0,9]/float(newdata[0,7])
    print(newdata)
    np.savetxt(target_file, newdata, fmt='%.6f',delimiter=' ', newline='\n', encoding=None)

pavg_list = [float(x)*0.001 for x in range(1,11)]
range_list=[0,0.1,0.2,0.3,0.6,0.9]
for p in pavg_list:
    for r in range_list:
        combine( "n200_range%.1f_p%.3f_serial2.data "%(r,p),"n200_range%.1f_p%.3f_serial2_sum.data "%(r,p))
        combine( "n200_range%.1f_p%.3f_parellel2.data "%(r,p),"n200_range%.1f_p%.3f_parallel2_sum.data "%(r,p))
