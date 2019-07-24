import time
import timeit
print "here"
s = """
x = []
for i in range(1000):
	x.append(i**2)"""

t = timeit.Timer(stmt=s)
print t.timeit(number=1)
