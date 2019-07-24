import pycuda.driver as cuda
import pycuda.autoinit
#cuda.init()

kernel_source = """
__global__ void
cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz)
{
}
"""
cuda_texture=cuda.SourceModule(kernel_source)
