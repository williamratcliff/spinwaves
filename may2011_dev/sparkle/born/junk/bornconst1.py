# This program is in the public domain
# Authors: Paul Kienzle, Christopher Metting
#03/23/2010

import Queue
import threading
from pycuda import gpuarray
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy
import time
cuda.init()
def cudoSMBA_form(cell,Q,lattice,beam):
    '''
    This module ties the cudo calculation of the SMBA to my class structures. Other tie-ins may
    easily be included here.
    '''
    form =(born(cell.unit,cell.value_list[0],cell.value_list[1],cell.value_list[2],
                Q.q_list[0],Q.q_list[1],Q.q_list[2]))
    return form

def main():
    #unit_size = (50,50,50)
    #feature_size = (5,5,5)
    #Q_size = (12,11,10)
    unit_size = (50,50,50)
    unit_metric = (1.0e5,1.0e5,2.5e3)
    feature_size = (25,25,41)
    #Q_size = (300,300,300)
    Q_size = (100,100,100)
    #Q_size = (50,50,50)

    #print "done!",form[:,-1,-1].real+1
    #print form
    unit = numpy.zeros(unit_size, 'd')
    unit[0:feature_size[0],0:feature_size[1],0:feature_size[2]] = 4.5e-6
    gpu=None
    #gpu=1  # Use this line if you only want one GPU to be in use
    #gpu=-1  # Use this line if you want to leave one GPU free
    
    qx = numpy.linspace(-0.0001,0.0001,Q_size[0])
    qy = numpy.linspace(-0.0001,0.0001,Q_size[1])
    qz = numpy.linspace(0.00002,0.2,Q_size[2])

    x = (unit_metric[0]/unit_size[0])*numpy.arange(unit_size[0])
    y = (unit_metric[1]/unit_size[1])*numpy.arange(unit_size[1])
    z = (unit_metric[2]/unit_size[2])*numpy.arange(unit_size[2])

    t0 = time.time()
    form32 = born(unit, x, y, z, qx, qy, qz, gpu=gpu, precision='float32')
    print "time",time.time()-t0
    #print "count nans 32",numpy.sum(numpy.isnan(form32))
    ysum32 = (numpy.sum(abs(form32)**2, axis=1))
    ysum32 += numpy.min(ysum32[(ysum32)>0])/2


    if 1:
        import pylab
        pylab.figure(0)
        pylab.imshow(numpy.log10(ysum32))
        pylab.colorbar()

    if 0:
        form64 = born(unit, x, y, z, qx, qy, qz, gpu=gpu, precision='float64')
        #print "count nans 64",numpy.sum(numpy.isnan(form64))
        ysum64 = (numpy.sum(abs(form64)**2, axis=1))
        ysum64 += numpy.min(ysum64[(ysum64)>0])/2
    
        print "max diff 64-32",abs(ysum64-ysum32).max()    

        pylab.figure(1)
        pylab.imshow(numpy.log10(ysum64))
        pylab.colorbar()
        pylab.figure(2)
        pylab.imshow((ysum64-ysum32)/ysum64)
        pylab.colorbar()
    
    pylab.show()
    

# Note: may want to try putting const arrays in texture memory rather
# than global memory; claims that card can do better memory caching
# hence maybe get better performance.

def loadkernelsrc(name, precision='float32', defines={}):
    import os
    srcfile = os.path.join(os.path.dirname(__file__),name)
    file = open(srcfile)
    src = file.read()
    file.close()
    defines = "\n".join(('#define %s %s'%(k,str(defines[k]))) 
                        for k in sorted(defines.keys()))
        
    if precision == 'float32':
        typedefs = '''
typedef float real;
typedef float2 cplx;
#define sin __sinf
#define cos __cosf
#define sincos __sincosf
        '''
    else:
        typedefs = '''
typedef double real;
typedef double2 cplx;
        '''
    src = defines+typedefs+src
    return src

#BORN_KERNEL32 = loadkernelsrc('kernelconst.c', precision='float32')
#BORN_KERNEL64 = loadkernelsrc('kernelconst.c', precision='float64')

def born(density, x, y, z, Qx, Qy, Qz, gpu=None, precision='float32'):
    density,x,y,z,Qx,Qy,Qz = [numpy.asarray(v,precision) for v in density,x,y,z,Qx,Qy,Qz]
    cplx = 'complex64' if precision=='float32' else 'complex128'
    if gpu is -1:
        numgpus = cuda.Device.count() - 1
        gpus = range(1,numgpus+1)
    elif gpu is not None:
        numgpus = 1
        gpus = [gpu]
    else:
        numgpus = cuda.Device.count()
        gpus = range(numgpus)
    size = [len(v) for v in Qx,Qy,Qz]
    result = numpy.empty(size,dtype=cplx)
    work_queue = Queue.Queue()
    for qxi,qx in enumerate(Qx): work_queue.put(qxi)
    

    threads = [BornThread(gpus[k], work_queue, result,
                      density, x, y, z, Qx, Qy, Qz)
           for k in range(numgpus)]

    
    for T in threads: T.start()
    for T in threads: T.join()
    return result

class BornThread(threading.Thread):
    def __init__(self, gpu, work_queue, result, 
                 density, x, y, z, Qx, Qy, Qz):
        threading.Thread.__init__(self)
        self.born_args = density, x, y, z, Qx, Qy, Qz, result
        self.work_queue = work_queue
        self.gpu = gpu
        self.precision = Qx.dtype
        self.defines = dict(XSIZE=len(x),YSIZE=len(y),ZSIZE=len(z),
                            QXSIZE=len(Qx),QYSIZE=len(Qy),QZSIZE=len(Qz))
        
    def run(self):
        self.dev = cuda.Device(self.gpu)
        self.ctx = self.dev.make_context()
        src = loadkernelsrc("kernelconst1.c",
                            precision=self.precision,
                            defines=self.defines)
        #print src
        self.cudamod = SourceModule(src)
        self.cudaBorn = self.cudamod.get_function("cudaBorn")
        self.kernel()
        self.ctx.pop()
        del self.ctx
        del self.dev


    def kernel(self):
        density, x, y, z, Qx, Qy, Qz, result = self.born_args
        nx, ny, nz, nqx, nqy, nqz = [numpy.int32(len(v)) 
                                     for v in x, y, z, Qx, Qy, Qz]
        freelist = []
        for v,t in [(density,"tdensity")]:
            cv = gpuarray.to_gpu(v.flatten())
            tv = self.cudamod.get_texref(t)
            cv.bind_to_texref_ext(tv)
            freelist.append(cv)
        for v,t in (x,"x"),(y,"y"),(z,"z"),(Qx,"Qx"),(Qy,"Qy"),(Qz,"Qz"):
            cv = self.cudamod.get_global(t)[0]
	    cuda.memcpy_htod(cv, v)
            freelist.append(cv)

        cframe = cuda.mem_alloc(result[0].nbytes)
        n = int(1*nqy*nqz)
        
        while True:
            try:
                qxi = numpy.int32(self.work_queue.get(block=False))
            except Queue.Empty:
                break
            
            #print "%d of %d on %d\n"%(qxi,nqx,self.gpu),
            
            self.cudaBorn(qxi, cframe, **cuda_partition(n))
            ## Delay fetching result until the kernel is complete
            cuda_sync()

            ## Fetch result back to the CPU
            cuda.memcpy_dtoh(result[qxi], cframe)

            #print "%d %s\n"%(qxi,ctemp.get())

def cuda_sync():
    """
    Waits for operation in the current context to complete.
    """
    #return # The following works in C++; don't know what pycuda is doing
    # Create an event with which to synchronize
    done = cuda.Event()
    # Schedule an event trigger on the GPU.
    done.record()
    # Block until the GPU executes the kernel; the sync call is a busy loop,
    # so put a timer in here to pause.
    while not done.query(): time.sleep(0.001)
    done.synchronize()
    # Event is deleted on completion of the function

def cuda_partition(n):
    max_gx,max_gy = 65535,65535
    blocksize = 32
    #max_gx,max_gy = 5,65536
    #blocksize = 3
    block = (blocksize,1,1)
    num_blocks = int((n+blocksize-1)/blocksize)
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


if __name__ == "__main__":
    main()
