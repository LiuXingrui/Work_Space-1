import os
os.system("rm wt_n200d.txt")
os.system("rm wt_n200s.txt")
a=range(1,7)
for n in a:
    command="./qBP toric_Hx_n10_cyc11 toric_Hz_n10_cyc11 %f 10 1000 20 wt_n200s.txt "%n

    os.system(command)
    command="./qBP toric_Hx_n10_cyc11 toric_Hz_n10_cyc11 %f -1 1000 20 wt_n200d.txt "%n
    os.system(command)
