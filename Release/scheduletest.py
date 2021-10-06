import os

complexity = [[10000,500],[100000,3500],[200000,5000],[500000,7000]]

#os.system('setx OMP_SCHEDULE "static,0"')

print("schedule difference")
os.system('echo ScheduleDifferences >> outSchedule.txt')
for c in complexity:
    N = c[0]
    N_m = c[1]
    for crs in [1,2,4,6,8,10,12,14,16,32]:
        os.system('ScheduleDifferences.exe ' + str(N) + ' ' + str(crs) + ' 1 >> outSchedule.txt')
    os.system('echo " " >> outSchedule.txt')

print("Done")
