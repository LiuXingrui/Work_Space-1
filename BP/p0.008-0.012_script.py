import os
os.system("rm p0.008-0.012_data.txt")
a=range(3,18)
for n in a:
    command4="./qBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.008 0.012 2000 20 p0.008-0.012_data.txt"%(n,n)
    os.system(command4)
