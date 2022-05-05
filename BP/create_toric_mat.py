import os

a=range(3,21)
for n in a:
   # subprocess.run(["./gen_cyc2", "testn%d_cyc11"%n,"%d"%n,"1 1"]) #subprocess will give a wrong file, don't know why
    command1="./gen_cyc n%d_cyc11 %d 1 1"%(n, n)
    os.system(command1)
    #print(command1)
    command2="./tran n%d_cyc11 n%d_cyc11T"%(n,n)
    #print(command2)
    os.system(command2)
    command3="./HPC  n%d_cyc11 n%d_cyc11T  toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11"%(n,n,n,n)
    os.system(command3)
   # print(command3)
   # command4="./qBP toric_Hx_n%d_cyc11 toric_Hz_n%d_cyc11 0.01 0.01 2000 20 data2.txt"%(n,n)
   # os.system(command4)
   # print(command4)
    

