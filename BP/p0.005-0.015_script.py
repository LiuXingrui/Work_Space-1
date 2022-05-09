import os
os.system("rm p0.005-0.015_data.txt")
a=range(3,18)
for n in a:
    command4="./qBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.005 0.015 2000 20 p0.005-0.015_data.txt"%(n,n)
    os.system(command4)
