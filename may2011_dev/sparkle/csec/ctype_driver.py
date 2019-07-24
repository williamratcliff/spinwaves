"""
I am calling my CUDA code from c and just packaging everything in a shared object.
Then I will use ctypes to call the driver function(s) from python.

To compile the shared object file, I used the command:

nvcc --shared -o _mcqueenylib.so McQueeny_Alg.cu --compiler-options '-fPIC'
"""

import os
import sys
import ctypes
from ctypes import c_float, c_int, c_long
import time
#import matplotlib.pyplot as plt
import numpy as np
import timeit


#I'm not really sure why, but a SegFault was created when python exited
#it had to do with the cleanup after calling the C code, and using
#pycuda's initialization fixes it for some reason.
import pycuda.autoinit

def loadLib():
    if sys.platform in ('darwin'):
        ext = '.dylib'
    elif sys.platform in ('win32','cygwin'):
        ext = '.pyd'
    else:
        ext = '.so'

    dllpath=os.path.join(os.path.dirname(__file__),'_test'+ext)#if executing .py files
    dllpath = "/home/tsarvey/csec/_mcqueenylib.so"
    print "path: ", dllpath
    mcqueenyDll = ctypes.CDLL(dllpath)
    return mcqueenyDll

mcQueenyDll = loadLib()

mcQueenyDll.test_func(5)


def func(i):
	mcQueenyDll.test_func(i)

times = []
x = range(0,10000,1)
s1 = "func("
s2 = ")"
for i in x:
	s =  s1 + str(i) + s2
	print "string: ", s
	#t = timeit.Timer(s, "from __main__ import func")
	#times.append(t.timeit(number=1))
	mcQueenyDll.cuda_test(i)

#print "times: ", times
#plt.plot(x,times)
#plt.show()
np.savetxt("times.txt", times)

