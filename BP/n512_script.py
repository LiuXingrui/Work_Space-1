import os
os.system("rm n512_data.txt")
a=range(1,4)
for n in a:
    command4="./qBP toric_Hx_n16_cyc11 toric_Hz_n16_cyc11 0.005 0.015 10 20 n392_data.txt "
    command5="./qBP toric_Hx_n13_cyc11 toric_Hz_n13_cyc11 0.005 0.015 10 20 n392_data.txt "
    os.system(command4)
   # os.system(command5)
