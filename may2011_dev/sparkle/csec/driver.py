import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy
import time
cuda.init()


#Taken from Paul's born.py:
# Note: may want to try putting const arrays in texture memory rather
# than global memory; claims that card can do better memory caching
# hence maybe get better performance.

def readfile(name):
    file = open(name)
    txt = file.read()
    file.close()
    return txt

#Taken from Paul's born.py:
def loadkernel(name, precision='float32', defines={}):
    import os
    src = readfile(os.path.join(os.path.dirname(__file__),name))
    defines = "\n".join(('#define %s %s'%(k,str(defines[k]))) 
                        for k in sorted(defines.keys()))
        
    if precision == 'float32':
        typedefs = '''
#define HAVE_CUDA
#include <cufloat.h>
#define sin __sinf
#define cos __cosf
        '''
    else:
        typedefs = '''
#define HAVE_CUDA
#include <cudouble.h>
        '''
    src = defines+typedefs+src

    return SourceModule(src, no_extern_c=True, include_dirs=[os.path.abspath(os.path.dirname(__file__))])


def main():
	kernel = loadkernel("McQueeny_alg.cu")
	test_func = kernel.get_function("test_func");
	test_func()


if __name__ == "__main__":
    main()