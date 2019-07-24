import ctypes

#I'm not really sure why, but a SegFault was created when python exitted
#it had to do with the cleanup after calling the C code, and using
#pycudas initialization fixes it for some reason.
import pycuda.autoinit


mcQueenyDll = ctypes.CDLL("/home/tsarvey/csec/ctypes_test/_test.so")
mcQueenyDll.test_func()
print 'done'