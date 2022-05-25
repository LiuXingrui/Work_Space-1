import os
os.system("rm wt_n18da.txt")
os.system("rm wt_n18sa.txt")
a=range(2,3)
for n in a:
    print("wt 2, debug=0, range=0")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 2 0 1000 100 wt_n18sa.txt 0"
    os.system(command)

    
    
    print("pavg 0.01, debug=0, range=0")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 0.01 0 1000 100 wt_n18sa.txt 0"
    os.system(command)

    print("pavg 0.01, debug=0, range=0.5")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 0.01 0.5 1000 100 wt_n18sa.txt 0"
    os.system(command)

    print("wt 2, debug=1, range=0.5")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 2 0.5 1000 100 wt_n18sa.txt 1"
    os.system(command)

    print("pavg 0.01, debug=2,parallel, range=0.5")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 0.01 0.5 1000 100 wt_n18sa.txt 2"
    os.system(command)

    print("wt 2, debug=2,parallel, range=0.5")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 2 0.5 1000 100 wt_n18sa.txt 2"
    os.system(command)


    print("wt 2, debug=7,para ana, range=0.5")
    command="./xqBP toric_Hx_n3_cyc11 toric_Hz_n3_cyc11 2 0.5 10 10 wt_n18sa.txt 7"
    os.system(command)
    
 
