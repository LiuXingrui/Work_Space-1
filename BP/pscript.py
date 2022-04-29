import os

a=range(3,13)
for n in a:
    command4="./xqBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.01 2000 20"%(n,n)
    os.system(command4)
   # print(command4)
    

