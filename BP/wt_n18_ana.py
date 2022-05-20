import os
os.system("rm wt_n18da.txt")
os.system("rm wt_n18sa.txt")
a=range(2,3)
for n in a:
    command="./aqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 %f 0 50 100 wt_n18sa.txt "%n

    os.system(command)
    command="./aqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 %f 0.2 50 100 wt_n18da.txt "%n
    os.system(command)
