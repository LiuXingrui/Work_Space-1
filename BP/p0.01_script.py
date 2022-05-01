import os
os.system("rm p0.01_data.txt")
a=range(3,10)
for n in a:
    command4="./qBP toric_Hx_n12_cyc11 toric_Hz_n12_cyc11 %f 100 20 n12_data_sp.txt"%(0.001*n)
    os.system(command4)
