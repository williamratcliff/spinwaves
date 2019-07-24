import pycuda.driver as cuda
import pycuda.autoinit
import numpy

realrow = numpy.array([1.0, 2.0, 3.0, 4.0, 5.0],
                dtype=numpy.float32).reshape(1,5)

mod_copy_texture=cuda.SourceModule("""
texture<float, 1> tex;
texture<float, 1> tex2;

__global__ void copy_texture_kernel(float * data) {
    int ty=threadIdx.y;
    //data[ty] = tex1D(tex, (float)(ty));
    data[ty] = tex1D(tex, (float)(ty)/2.0f);
}
""")

copy_texture_func = mod_copy_texture.get_function("copy_texture_kernel")
texref = mod_copy_texture.get_texref("tex")
tex2ref = mod_copy_texture.get_texref("tex2")

cuda.matrix_to_texref(realrow, texref, order="C")

texref.set_flags(cuda.TRSF_NORMALIZED_COORDINATES)
texref.set_filter_mode(cuda.filter_mode.LINEAR)

gpu_output = numpy.zeros_like(realrow)
copy_texture_func(cuda.Out(gpu_output), block=(1,1,1), texrefs=[texref])

print "Output:"
print gpu_output

