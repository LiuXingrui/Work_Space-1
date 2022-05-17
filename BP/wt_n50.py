import os
os.system("rm wt_n50d.txt")
os.system("rm wt_n50s.txt")
a=range(1,6)
for n in a:
    command="./qBP toric_Hx_n5_cyc11 toric_Hz_n5_cyc11 %f 10 1000 20 wt_n50s.txt "%n

    os.system(command)
    command="./qBP toric_Hx_n5_cyc11 toric_Hz_n5_cyc11 %f -1 1000 20 wt_n50d.txt "%n
    os.system(command)
