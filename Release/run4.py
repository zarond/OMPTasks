import os

#complexity = [[10000,500],[500000,3500],[1000000,5000],[2000000,7000]]

#os.system('set OMP_SCHEDULE=static,0')

print("Task4")
os.system('echo Task4 >> out4.txt')

for crs in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,25,32,36]:
    os.system('Task4.exe ' + str(crs) + ' ' + str(crs) + ' 1 >> out4.txt')


print("Done")
