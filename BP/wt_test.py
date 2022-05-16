import os
os.system("rm test.txt")
os.system("rm test2.txt")
a=range(1,6)
for n in a:
    command4="./qBP toric_Hx_n6_cyc11 toric_Hz_n6_cyc11 %f 10 300 20 test.txt "%n
    #command5="./qBP toric_Hx_n13_cyc11 toric_Hz_n13_cyc11 0.005 0.015 10 20 n392_data.txt "
    os.system(command4)
    command4="./qBP toric_Hx_n6_cyc11 toric_Hz_n6_cyc11 %f -1 300 20 test2.txt "%n
    os.system(command4)
