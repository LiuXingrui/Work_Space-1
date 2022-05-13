import os
os.system("rm n18_data.txt")
a=range(1,25)
for n in a:
    command4="./qBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 %f %f 2000 20 n18_data.txt "%(0.005*n,0.005*n)
    os.system(command4)
