# This program is in the public domain
# Authors: Paul Kienzle, Christopher Metting
#03/25/2010

import math
import time
from pycuda import gpuarray
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy


def cuda_sync():
    """
    Waits for operation in the current context to complete.
    """
    #return # The following works in C++; don't know what pycuda is doing
    # Create an event with which to synchronize
    done = cuda.Event()
    # Schedule an event trigger on the GPU.
    done.record()
    # Block until the GPU executes the kernel.
    done.synchronize()
    # Clean up the event; I don't think they can be reused.
    del done

def cuda_partition(n, blocksize=32):
    max_gx,max_gy = 65535,65535
    #max_gx,max_gy = 5,65536
    #blocksize = 3
    block = (blocksize,1,1)
    num_blocks = int((nx+blocksize-1)/blocksize)
    if num_blocks < max_gx:
        grid = (num_blocks,1)
    else:
        gx = max_gx
        gy = (num_blocks + max_gx - 1) / max_gx
        if gy >= max_gy: raise ValueError("vector is too large")
        grid = (gx,gy)
    #print "block",block,"grid",grid
    #print "waste",block[0]*block[1]*block[2]*grid[0]*grid[1] - n
    return dict(block=block,grid=grid)




kernel = """
    /*
    Full calculation with 3-D block and 2-D grid.
    
    // get my block within a grid
    const int myblock=blockIdx.x+blockIdx.y*gridDim.x;
    // how big is each block within a grid
    const int blocksize=blockDim.x*blockDim.y*blockDim.z;
    // get thread within a block
    const int subthread=threadIdx.z*(blockDim.x*blockDim.y)+threadIdx.y*blockDim.x+threadIdx.x;
    // find my thread
    const int thread=myblock*blocksize+subthread;
    */


__global__ void
f64(int nx,double x[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= nx ) return;    

    x[idx] = idx; // sin(x[idx]);
}

__global__ void
f32(int nx,float x[],int ny,float y[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= nx ) return;

    const float x0 = x[idx];
    float result = 0;
//# pragma unroll 4
    for (int j=0; j < 1000; j++)  {
        for (int k=0; k < ny; k++) {
            result +=  x0*__sinf(y[k]);
        }
    }
    x[idx] = result;
}
"""


cuda.init()
dev = cuda.Device(1)
ctx = dev.make_context()
#import pycuda.autoinit

cudamod = SourceModule(kernel)

nx = 100000
ny = 100
if 1:
    t0 = time.clock()
    cudaf32 = cudamod.get_function("f32")
    x32 = numpy.arange(nx,dtype = 'float32')
    z32 = numpy.arange(ny,dtype = 'float32')
    cudax32 = gpuarray.to_gpu(x32)
    cudaz32 = gpuarray.to_gpu(z32)
    cudaf32(numpy.int32(nx), cudax32, 
            numpy.int32(ny), cudaz32,
            shared=32*4, **cuda_partition(nx))
    cuda_sync()
    y32 = cudax32.get()
    print "Time for 32 bit",time.clock() - t0

if 0:
    t0 = time.clock()
    cudaf64 = cudamod.get_function("f64")
    x64 = numpy.arange(nx,dtype = 'float64')
    cudax64 = gpuarray.to_gpu(x64)
    cudaf64(numpy.int32(nx),cudax64,
            **cuda_partition(nx))
    cuda_sync()
    y64 = cudax64.get()
    print "Time for 64 bit",time.clock() - t0

#print x32.dtype
#print y32.dtype

#print x64.dtype
#print y64.dtype


#print 'cuda32',y32
#print 'python', numpy.sin(x32)

#print 'cuda64',y64
#print 'python', numpy.sin(x64)

#print 'diffs', numpy.linalg.norm(numpy.sin(x32)-y32)
#print 'diffs', numpy.linalg.norm(numpy.sin(x64)-y64)
#print 'diffs', numpy.linalg.norm(numpy.sin(x32)-y64)
#print 'diffs', numpy.linalg.norm(numpy.sin(x64)-y32)
#print 'diffs', numpy.linalg.norm(y64-y32)
ctx.pop()
del ctx, dev
