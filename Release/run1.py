import os

complexity = [10000,50000,100000,500000,700000,1000000,1250000,1500000,1750000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000]

#os.system('set OMP_SCHEDULE=static,0')

print("Task1")
os.system('echo Task1 >> out1.txt')
for c in complexity:
    N = c
    for crs in [1,2,4,6,8,10,12,14,16,32]:
        os.system('Task1.exe ' + str(N) + ' ' + str(crs) + ' 1 >> out1.txt')
    os.system('echo " " >> out1.txt')



print("Done")
