texture <real, 3> tdensity;
__constant__ real Qx[QXSIZE], Qy[QYSIZE], Qz[QZSIZE], x[XSIZE], y[YSIZE], z[ZSIZE];
__global__ void
cudaBorn(int qxi, cplx result[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= QXSIZE*QYSIZE*QZSIZE) return;
    const int qzi = idx%QZSIZE;
    const int qyi = (idx/QZSIZE)%QYSIZE;

#if 1
    real Re = 0;
    real Im = 0;
    for (int xi=0; xi < XSIZE; xi++) {
        for (int yi=0; yi < YSIZE; yi++) {
#if 0 // No unrolling: 36 s (940 Gflops)
            for (int zi=0; zi < ZSIZE; zi++) {
                const real QdotR = x[xi]*Qx[qxi]+y[yi]*Qy[qyi]+z[zi]*Qz[qzi];
                const real rho = tex3D(tdensity,zi,yi,xi);
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
            }
#else  // Loop unrolling: 29 s (1160 Gflops) 630 W
            const real xyQxy = x[xi]*Qx[qxi] + y[yi]*Qy[qyi];
            int zi;
            for (zi=0; zi < ZSIZE-(ZSIZE%4); zi+=4) {
                const real QdotR0 = xyQxy+z[zi]*Qz[qzi];
                const real QdotR1 = xyQxy+z[zi+1]*Qz[qzi];
                const real QdotR2 = xyQxy+z[zi+2]*Qz[qzi];
                const real QdotR3 = xyQxy+z[zi+3]*Qz[qzi];
                const real rho0 = tex3D(tdensity,zi,yi,xi);
                const real rho1 = tex3D(tdensity,zi+1,yi,xi);
                const real rho2 = tex3D(tdensity,zi+2,yi,xi);
                const real rho3 = tex3D(tdensity,zi+3,yi,xi);
                Re += rho0*cos(QdotR0) + rho1*cos(QdotR1) + rho2*cos(QdotR2) + rho3*cos(QdotR3);
                Im += rho0*sin(QdotR0) + rho1*sin(QdotR1) + rho2*sin(QdotR2) + rho3*sin(QdotR3);
            }
            for (; zi < ZSIZE; zi++) {
                const real QdotR = xyQxy+z[zi]*Qz[qzi];
                const real rho = tex3D(tdensity,zi,yi,xi);
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
            }
#endif
        }
    }
    result[idx].x = Re;
    result[idx].y = Im;
#else
    result[idx].x = result[idx].y = 0;
    if (qxi < XSIZE && qyi < YSIZE && qzi < ZSIZE) result[idx].x = tex3D(tdensity,qxi,qyi,qzi);
#endif
}
