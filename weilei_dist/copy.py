import os

a=range(3,18)
for n in a:
  
    command1="cp ../S_BP/toric_Hx_n%d_cyc11 toric_Hx_n%d_cyc11"%(n,n)
    os.system(command1)
    command1="cp ../S_BP/toric_Hz_n%d_cyc11 toric_Hz_n%d_cyc11"%(n,n)
    os.system(command1)
