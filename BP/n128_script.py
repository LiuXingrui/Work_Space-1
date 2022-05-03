import os
os.system("rm n128_data.txt")
a=range(1,15)
for n in a:
    command4="./qBP toric_Hx_n8_cyc11 toric_Hz_n8_cyc11 %f %f 2000 20 n128_data.txt "%(0.004*n,0.006*n)
    os.system(command4)
