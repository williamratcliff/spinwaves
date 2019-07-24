#import pyximport; pyximport.install()
import hello


import time

init_time = time.clock()
j = 0
for i in range(0,100000):
    j+=i

print 'python'
print j
print time.clock()-init_time

