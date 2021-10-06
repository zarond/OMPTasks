import os

complexity = [[10000,500],[500000,3500],[1000000,5000],[2000000,7000]]

#os.system('setx OMP_SCHEDULE "static,0"')

print("Task5")
os.system('echo Task5 >> out.txt')
#os.system('echo Task5 static >> out.txt')
#os.system('set OMP_SCHEDULE=static,0')
for c in complexity:
    N = c[0]
    N_m = c[1]
    for crs in [1,2,4,6,8,10,12,14,16,32]:
        os.system('Task5.exe ' + str(N_m) + ' ' + str(crs) + ' 1 >> out.txt')
    os.system('echo " " >> out.txt')

print("Done")
