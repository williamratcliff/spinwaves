extern "C" void test_func()
{
    int *dev_a;

    // allocate the memory on the GPU (array of 2 ints)
    cudaMalloc( (void**)&dev_a, 2 * sizeof(int) );
    cudaFree( dev_a );   
}
