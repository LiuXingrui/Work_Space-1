import os
os.system("rm wt_n18d.txt")
os.system("rm wt_n18s.txt")
a=range(1,4)
for n in a:
    command="./qBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 %f 10 1000 20 wt_n18s.txt "%n

    os.system(command)
    command="./qBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 %f -1 1000 20 wt_n18d.txt "%n
    os.system(command)
