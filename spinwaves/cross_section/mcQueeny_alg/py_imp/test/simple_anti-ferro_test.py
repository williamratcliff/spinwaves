""" I am testing McQueeny's aglorithm on the simple anti-ferromagnetic case.
I need to check which way the indices on the eigenvectors go.  I already
calculated what the eigenvectors and eigenvalues should be in my notebook.
This produces a plot of S vs. q, without any gamma, r_0, k'/k, or deltas"""

import numpy as np
from numpy import sqrt, exp, sin
from matplotlib import pyplot as plt
import math




def dbl_equal(d1, d2, maxAbsoluteError = 1.0e-10, maxRelativeError = 1.0e-5):
    """A good discussion of floating point equality comparison can be found here:
    http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    
    I adapted this form one of his C AlmostEqual functions."""
    if d1 == d2:
        return True
    if math.isnan(d1.real) or math.isnan(d2.real) or math.isnan(d1.imag) or math.isnan(d2.imag):#In this case I want this
        print d1, " equal to ", d2
        return True
    if abs(d1-d2) < maxAbsoluteError:
        return True
    if abs(d2) > abs(d1):
        relativeError = abs((d1 - d2) / d2)
    else:
        relativeError = abs((d1 - d2) / d1)
    if relativeError <= maxRelativeError:
        return True;
    print d1, " not equal to ", d2, "    (relative error = ", relativeError, ")"
    return False






i = 1j
#This is only valid between 0 and 2pi
q_range = np.arange(0,2*np.pi,0.01)
S_q = [] 

for q in q_range:

    M = np.array([
    [2, -1-exp(-i*q)],
    [1+exp(i*q), -2]])
    
    V1 = np.array([
    [( 2-sqrt(-exp(-i*q)*(-1+exp(i*q))**2) )/(1+exp(i*q))],
    [1]])
    
    V2 = np.array([
    [( 2+sqrt(-exp(-i*q)*(-1+exp(i*q))**2) )/(1+exp(i*q))],
    [1]])

    C1 = sqrt( -(1-sin(q/2)**2) / ( 2 * (sin(q/2)**2 - sin(q/2)) ) )
    
    C2 = sqrt( (1-sin(q/2)**2) / ( 2 * (sin(q/2)**2 + sin(q/2)) ) )
    
    T1 = C1*V1
    T2 = C2*V2
    
    #now let's check to make sure that these are properly normalized
    #and that they are indeed eigenvectors of M
    
    l1 = -2*sin(q/2)
    l2 = 2*sin(q/2)
    
    assert dbl_equal(l1*T1[0], np.dot(M,T1)[0])
    assert dbl_equal(l1*T1[1], np.dot(M,T1)[1])
    
    assert dbl_equal(l2*T2[0], np.dot(M,T2)[0])
    assert dbl_equal(l2*T2[1], np.dot(M,T2)[1])
    
    assert dbl_equal(abs(T1[0])**2 - abs(T1[1])**2, -1.0)
    assert dbl_equal(abs(T2[0])**2 - abs(T2[1])**2, 1.0)
    
    #All assertions are good
    #correct
    S_q.append( (abs( T1[0]*exp(-i*q/4) - T1[1]*exp(-i*3*q/4) )**2 + abs( T2[0]*exp(-i*q/4) - T2[1]*exp(-i*3*q/4) )**2)[0] )
    #false
    #S_q.append( (abs( T1[0]*exp(-i*q/4) - T2[0]*exp(-i*3*q/4) )**2 + abs( T1[1]*exp(-i*q/4) - T2[1]*exp(-i*3*q/4) )**2)[0] )
    #print "\nq: ", q
    #print "T1: ", T1, "   w1: ", l1 
    #print "T2: ", T2, "   w2: ", l2
    #print "S_q: ", S_q[-1]
    print q, "\t", S_q[-1]
    

plt.plot(q_range[:-100], S_q[:-100])
plt.show()
