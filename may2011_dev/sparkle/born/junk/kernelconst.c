texture <real, 1> tdensity;
__constant__ real Qx[QXSIZE], Qy[QYSIZE], Qz[QZSIZE], x[XSIZE], y[YSIZE], z[ZSIZE];
__global__ void
cudaBorn(int qxi, cplx result[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= QXSIZE*QYSIZE*QZSIZE) return;
    const int qzi = idx%QZSIZE;
    const int qyi = (idx/QZSIZE)%QYSIZE;

    real Re = 0;
    real Im = 0;
    int densityidx = 0;
    for (int xi=0; xi < XSIZE; xi++) {
        for (int yi=0; yi < YSIZE; yi++) {
#if 0 // No unrolling: 36 s (940 Gflops)
            for (int zi=0; zi < nz; zi++) {
                const real QdotR = x[xi]*Qx[qxi]+y[yi]*Qy[qyi]+z[zi]*Qz[qzi];
                const real rho = tex1Dfetch(tdensity,densityidx);
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
                densityidx++;
            }
#else  // Loop unrolling: 29 s (1160 Gflops) 630 W
            const real xyQxy = x[xi]*Qx[qxi] + y[yi]*Qy[qyi];
            int zi;
            for (zi=0; zi < ZSIZE-(ZSIZE%4); zi+=4) {
                const real QdotR0 = xyQxy+z[zi]*Qz[qzi];
                const real QdotR1 = xyQxy+z[zi+1]*Qz[qzi];
                const real QdotR2 = xyQxy+z[zi+2]*Qz[qzi];
                const real QdotR3 = xyQxy+z[zi+3]*Qz[qzi];
                const real rho0 = tex1Dfetch(tdensity,densityidx);
                const real rho1 = tex1Dfetch(tdensity,densityidx+1);
                const real rho2 = tex1Dfetch(tdensity,densityidx+2);
                const real rho3 = tex1Dfetch(tdensity,densityidx+3);
                Re += rho0*cos(QdotR0) + rho1*cos(QdotR1) + rho2*cos(QdotR2) + rho3*cos(QdotR3);
                Im += rho0*sin(QdotR0) + rho1*sin(QdotR1) + rho2*sin(QdotR2) + rho3*sin(QdotR3);
                densityidx+=4;
            }
# if 1
            for (; zi < ZSIZE; zi++) {
                const real QdotR = xyQxy+z[zi]*Qz[qzi];
                const real rho = tex1Dfetch(tdensity,densityidx);
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
                densityidx++;
            }
# else  
# error "This is faster but for some reason doesn't produce the correct answer." 
// 26 s (1300 Gflops) 660 W
#  if ZSIZE%4 == 1
            const real QdotR0 = xyQxy+z[zi]*Qz[qzi];
            const real rho0 = tex1Dfetch(tdensity,densityidx);
            Re += rho0*cos(QdotR0);
            Im += rho0*sin(QdotR0);
#  endif
#  if ZSIZE%4 == 2
            const real QdotR0 = xyQxy+z[zi]*Qz[qzi];
            const real QdotR1 = xyQxy+z[zi+1]*Qz[qzi];
            const real rho0 = tex1Dfetch(tdensity,densityidx);
            const real rho1 = tex1Dfetch(tdensity,densityidx+1);
            Re += rho0*cos(QdotR0) + rho1*cos(QdotR1);
            Im += rho0*sin(QdotR0) + rho1*sin(QdotR1);
#  endif
#  if ZSIZE%4 == 3 
            const real QdotR0 = xyQxy+z[zi]*Qz[qzi];
            const real QdotR1 = xyQxy+z[zi+1]*Qz[qzi];
            const real QdotR2 = xyQxy+z[zi+2]*Qz[qzi];
            const real rho0 = tex1Dfetch(tdensity,densityidx);
            const real rho1 = tex1Dfetch(tdensity,densityidx+1);
            const real rho2 = tex1Dfetch(tdensity,densityidx+2);
            Re += rho0*cos(QdotR0) + rho1*cos(QdotR1) + rho2*cos(QdotR2);
            Im += rho0*sin(QdotR0) + rho1*sin(QdotR1) + rho2*sin(QdotR2);
#  endif
# endif
#endif
        }
    }
    result[idx].x = Re;
    result[idx].y = Im;
}
