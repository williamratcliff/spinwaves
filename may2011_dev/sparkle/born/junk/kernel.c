
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
    real Re = 0;
    real Im = 0;
    int densityidx = 0;
    for (int xi=0; xi < nx; xi++) {
        for (int yi=0; yi < ny; yi++) {
# pragma unroll 4
            for (int zi=0; zi < nz; zi++) {
                const real QdotR = Qx[qxi]*x[xi] + Qy[qyi]*y[yi] + Qz[qzi]*z[zi];
                Re += density[densityidx]*cos(QdotR);
                Im += density[densityidx]*sin(QdotR);
                densityidx++;
            }
        }
    }
    result[idx].x = Re;
    result[idx].y = Im;
}
