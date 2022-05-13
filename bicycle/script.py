import os
#os.system("rm ../../../Documents/Work_Space-1/S_BP/bicycle_Hx1") 
#os.system("rm ../../../Documents/Work_Space-1/S_BP/bicycle_Hz1")

a=range(1,10)
for i in a:
    t2=10*i
    command4="./bike.out  ../../../Documents/Work_Space-1/S_BP/bicycle_Hx%d  ../../../Documents/Work_Space-1/S_BP/bicycle_Hz%d %d 3 3 0 1 3 0 4 17"%(i,i,t2*2+1)
    os.system(command4)
