#define NWARP 1
#define WARP 32
#define BLK (WARP*NWARP)

__global__ void
cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz,
		const real density[],
		const real x[], const real y[], const real z[],
		const real Qx[], const real Qy[], const real Qz[],
		int qxi, cplx result[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx > nqx*nqy*nqz) return;
    const int qzi = idx%nqz;
    const int qyi = (idx/nqz)%nqy;
    const real qxs = Qx[qxi], qys = Qy[qyi], qzs = Qz[qzi];

    real Re = 1;
    real Im = 1;
    int size = nx*ny*nz;
    int di;
    __shared__ real xs[BLK], ys[BLK], zs[BLK], ds[BLK];
    int tail = size%BLK;
    size -= tail;
    // Use shared memory to cach x,y,z and density
    for (di = 0; di < size; di += BLK) {
        // Pull section of x,y,z,density into shared memory
        for (int i=0; i < NWARP; i++) {
            const int k = i + NWARP*threadIdx.x;
            xs[k] = x[(di+k)/(ny*nz)];
            ys[k] = y[((di+k)/nz)%ny];
            zs[k] = z[(di+k)%nz];
            ds[k] = density[di+k];
        }
        syncthreads();
        // Compute effects of x,y,z,density section
        for (int k=0; k < BLK; k++) {
            const real QdotR = xs[k]*qxs + ys[k]*qys + zs[k]*qzs;
            Re += ds[k]*cos(QdotR);
            Im += ds[k]*sin(QdotR);
        }
        syncthreads();
    }
    // Nothing fancy for the tail.
    while (tail) {
        const real QdotR = x[di/(ny*nz)]*qxs + y[(di/nz)%ny]*qys + z[di%nz]*qzs;
        Re += density[di]*cos(QdotR);
        Im += density[di]*sin(QdotR);
        di++; tail--;
    }

    result[idx].x = Re;
    result[idx].y = Im;
}
