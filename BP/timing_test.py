import os

a=range(1,2)
for n in a:
    command="./xqBP toric_Hx_n17_cyc11 toric_Hz_n17_cyc11 0.01 0.2 10000 100 test.data 1 0 1 0.05 0.2"
    os.system(command)

    command="./hqBP toric_Hx_n17_cyc11 toric_Hz_n17_cyc11 0.01 0.2 10000 100 test.data 1 0 1 0.05 0.2"
    os.system(command)

