import os
os.system("rm n128sp_data.txt")
a=range(1,15)
for n in a:
    command4="./qBP toric_Hx_n8_cyc11 toric_Hz_n8_cyc11 %f %f 2000 20 n128sp_data.txt"%(0.005*n,0.005*n)
    os.system(command4)
