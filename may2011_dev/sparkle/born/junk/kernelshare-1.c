#define WARP 32

__global__ void
cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz,
		const real density[],
		const real x[], const real y[], const real z[],
		const real Qx[], const real Qy[], const real Qz[],
		int qxi, cplx result[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    const int qzi = idx%nqz;
    const int qyi = (idx/nqz)%nqy;
    if (qyi >= nqy || qzi >= nqz) return;

    int densityidx = 0;
    real Re = 0;
    real Im = 0;
    __shared__ real zshared[WARP],dshared[WARP];
    for (int xi=0; xi < nx; xi++) {
        for (int yi=0; yi < ny; yi++) {
            for (int zi=0; zi < nz; zi+=WARP) {
                if (zi+threadIdx.x < nz) {
                    zshared[threadIdx.x] = z[zi+threadIdx.x];
                    dshared[threadIdx.x] = density[densityidx+threadIdx.x];
    	            const real xyQxy = x[xi]*Qx[qxi] + y[yi]*Qy[qyi];
                    const real qz = Qz[qzi];
                    const int end = (zi+WARP < nz) ? WARP : nz-zi;
                    syncthreads();
                    for (int w=0; w < nz; w++) {
                        const real QdotR = xqx + yqy + qz*zshared[w];
                        real cx,sx;
                        sincos(QdotR,&sx,&cx);
                        Re += dshared[w]*cx;
                        Im += dshared[w]*sx;
                    }
                }
                densityidx+=WARP;
            }
        }
    }
    result[idx].x = Re;
    result[idx].y = Im;
}
