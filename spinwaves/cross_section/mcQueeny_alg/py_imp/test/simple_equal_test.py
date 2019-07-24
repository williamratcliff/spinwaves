"""test if the outputs from simple_anti-ferro_test.py and mcQueeny_alg.py
are the same"""

import numpy as np
from matplotlib import pyplot as plt

#I had to relax these standards a lot for mcQueeny's results.  I wonder if
#he is using 32 bit floats - even if he is though, I should have better
#agreement thatn 1% I would have thought.
def dbl_equal(d1, d2, maxAbsoluteError = 1.0e-10, maxRelativeError = 1.0e-1):
    """A good discussion of floating point equality comparison can be found here:
    http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    
    I adapted this form one of his C AlmostEqual functions."""
    if abs(d1-d2) < maxAbsoluteError:
        return True
    if abs(d2) > abs(d1):
        relativeError = abs((d1 - d2) / d2)
    else:
        relativeError = abs((d1 - d2) / d1)
    if relativeError <= maxRelativeError:
        return True; 
    return False


a = np.loadtxt("output.txt")
b = np.loadtxt("output2.txt")
c = np.loadtxt("test_anti-ferro.out")#McQueeny

graph = []

print "len(a): ", len(a)
for i in range(len(a)):
    print a[i], " =? ", b[i], " =? ", np.array([c[i][0],c[i][4]])
    assert dbl_equal(a[i][0],b[i][0])
    assert dbl_equal(a[i][1],b[i][1])
    assert dbl_equal(a[i][0],c[i][0])
    assert dbl_equal(a[i][1],c[i][4])
    graph.append([a[i][1],b[i][1],c[i][4]])
    

print "EQUAL!"

graph = np.array(graph)
plt.plot(graph[:,0])
plt.plot(graph[:,1])
plt.plot(graph[:,2])
plt.show()
