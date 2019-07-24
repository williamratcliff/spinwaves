__constant__ real Qx[QXSIZE], Qy[QYSIZE], Qz[QZSIZE], x[XSIZE], y[YSIZE], z[ZSIZE];
__global__ void
cudaBorn(const real density[XSIZE][YSIZE][ZSIZE], int qxi, cplx result[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= QXSIZE*QYSIZE*QZSIZE) return;
    const int qzi = idx%QZSIZE;
    const int qyi = (idx/QZSIZE)%QYSIZE;

    real Re = 0;
    real Im = 0;
    for (int xi=0; xi < XSIZE; xi++) {
        for (int yi=0; yi < YSIZE; yi++) {
#if 0
# pragma unroll 4
            for (int zi=0; zi < ZSIZE; zi++) {
                const real QdotR = x[xi]*Qx[qxi]+y[yi]*Qy[qyi]+z[zi]*Qz[qzi];
                const real rho = density[xi][yi][zi];
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
            }
#else
            const real xyQxy = x[xi]*Qx[qxi] + y[yi]*Qy[qyi];
            int zi;
            for (zi=0; zi < ZSIZE-(ZSIZE%4); zi+=4) {
                const real QdotR0 = xyQxy+z[zi]*Qz[qzi];
                const real QdotR1 = xyQxy+z[zi+1]*Qz[qzi];
                const real QdotR2 = xyQxy+z[zi+2]*Qz[qzi];
                const real QdotR3 = xyQxy+z[zi+3]*Qz[qzi];
                const real rho0 = density[xi][yi][zi];
                const real rho1 = density[xi][yi][zi+1];
                const real rho2 = density[xi][yi][zi+2];
                const real rho3 = density[xi][yi][zi+3];
                Re += rho0*cos(QdotR0) + rho1*cos(QdotR1) + rho2*cos(QdotR2) + rho3*cos(QdotR3);
                Im += rho0*sin(QdotR0) + rho1*sin(QdotR1) + rho2*sin(QdotR2) + rho3*sin(QdotR3);
            }
            for (; zi < ZSIZE; zi++) {
                const real QdotR = xyQxy+z[zi]*Qz[qzi];
                const real rho = density[xi][yi][zi];
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
            }
#endif
        }
    }
    result[idx].x = Re;
    result[idx].y = Im;
}
