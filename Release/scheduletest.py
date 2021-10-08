import os

complexity = [[10000,500],[100000,3500],[200000,5000],[500000,7000]]
#chunks = [1,2,4,6,8,10,12,14,16,32]
chunks = [64,128,256,512,1024]

##print("schedule difference static")
##os.system('echo ScheduleDifferences static >> outSchedule.txt')
##for c in complexity:
##    N = c[0]
##    N_m = c[1]
##    for crs in [1,2,4,6,8,10,12,14,16,32]:
##        os.system('set OMP_SCHEDULE=static,0' + ' & ScheduleDifferences.exe ' + str(N) + ' ' + str(crs) + ' 1 >> outSchedule.txt')
##    os.system('echo " " >> outSchedule.txt')

print("schedule difference dynamic")

for ch in chunks:
    os.system('echo ScheduleDifferences dynamic,'+str(ch)+' >> outSchedule.txt')
    for c in complexity:
        N = c[0]
        N_m = c[1]
        for crs in [1,2,4,6,8,10,12,14,16,32]:
            os.system('set OMP_SCHEDULE=dynamic,'+str(ch) + ' & ScheduleDifferences.exe ' + str(N) + ' ' + str(crs) + ' 1 >> outSchedule.txt')
        os.system('echo " " >> outSchedule.txt')

print("schedule difference guided")

for ch in chunks:
    os.system('echo ScheduleDifferences guided,'+str(ch)+' >> outSchedule.txt')
    for c in complexity:
        N = c[0]
        N_m = c[1]
        for crs in [1,2,4,6,8,10,12,14,16,32]:
            os.system('set OMP_SCHEDULE=guided,'+str(ch) + ' & ScheduleDifferences.exe ' + str(N) + ' ' + str(crs) + ' 1 >> outSchedule.txt')
        os.system('echo " " >> outSchedule.txt')

print("Done")
