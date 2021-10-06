import os

complexity = [[10000,500],[500000,3500],[1000000,5000],[2000000,7000]]

#os.system('set OMP_SCHEDULE=static,0')

print("Task1")
os.system('echo Task1 >> out.txt')
for c in complexity:
    N = c[0]
    N_m = c[1]
    for crs in [1,2,4,6,8,10,12,14,16,32]:
        os.system('Task1.exe ' + str(N) + ' ' + str(crs) + ' 1 >> out.txt')
    os.system('echo " " >> out.txt')

print("Task2")
os.system('echo Task2 >> out.txt')
for c in complexity:
    N = c[0]
    N_m = c[1]
    for crs in [1,2,4,6,8,10,12,14,16,32]:
        os.system('Task2.exe ' + str(N) + ' ' + str(crs) + ' 1 >> out.txt')
    os.system('echo " " >> out.txt')

print("Task3")
os.system('echo Task3 >> out.txt')
for c in complexity:
    N = c[0]
    N_m = c[1]
    for crs in [1,2,4,6,8,10,12,14,16,32]:
        os.system('Task3.exe ' + str(N) + ' ' + str(crs) + ' 1 >> out.txt')
    os.system('echo " " >> out.txt')

print("Task4")
os.system('echo Task4 >> out.txt')
for c in complexity:
    N = c[0]
    N_m = c[1]
    for crs in [1,2,3,4,5,6,8,9,10,12,14,16,25,32,36]:
        os.system('Task4.exe ' + str(N_m) + ' ' + str(crs) + ' 1 >> out.txt')
    os.system('echo " " >> out.txt')

##print("Task5 static")
##os.system('echo Task5 static >> out.txt')
##os.system('set OMP_SCHEDULE=static,0')
##for c in complexity:
##    N = c[0]
##    N_m = c[1]
##    for crs in [1,2,4,6,8,10,12,14,16,32]:
##        os.system('Task5.exe ' + str(N_m) + ' ' + str(crs) + ' 1 >> out.txt')
##    os.system('echo " " >> out.txt')
##
##print("Task5 dynamic")
##os.system('echo Task5 dynamic >> out.txt')
##os.system('set OMP_SCHEDULE=dynamic,4')
##for c in complexity:
##    N = c[0]
##    N_m = c[1]
##    for crs in [1,2,4,6,8,10,12,14,16,32]:
##        os.system('Task5.exe ' + str(N_m) + ' ' + str(crs) + ' 1 >> out.txt')
##    os.system('echo " " >> out.txt')
##
##print("Task5 guided")
##os.system('echo Task5 guided >> out.txt')
##os.system('set OMP_SCHEDULE=guided,4')
##for c in complexity:
##    N = c[0]
##    N_m = c[1]
##    for crs in [1,2,4,6,8,10,12,14,16,32]:
##        os.system('Task5.exe ' + str(N_m) + ' ' + str(crs) + ' 1 >> out.txt')
##    os.system('echo " " >> out.txt')

##print("Task7")
##os.system('set OMP_SCHEDULE=static,0')
##os.system('echo Task7 >> out.txt')
##for c in complexity:
##    N = c[0]
##    N_m = c[1]
##    for crs in [1,2,4,6,8,10,12,14,16,32]:
##        os.system('Task7.exe ' + str(N) + ' ' + str(crs) + ' 1 >> out.txt')
##    os.system('echo " " >> out.txt')

print("Done")
