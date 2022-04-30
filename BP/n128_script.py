import os

a=range(1,21)
for n in a:
    command4="./qBP toric_Hx_n8_cyc11 toric_Hz_n8_cyc11 %f 2000 20 n8_data.txt"%(0.005*n)
    os.system(command4)
