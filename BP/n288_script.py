import os
os.system("rm n288_data.txt")
a=range(1,11)
for n in a:
    command4="./qBP toric_Hx_n12_cyc11 toric_Hz_n12_cyc11 %f 1000 20 n288_data.txt"%(0.005*n)
    os.system(command4)
