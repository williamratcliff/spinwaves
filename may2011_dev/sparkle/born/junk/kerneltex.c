texture <float, 1> tx, ty, tz, tQx, tQy, tQz, tdensity;
__global__ void
cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz, int qxi, cplx result[])
{
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= nqx*nqy*nqz) return;
    const int qzi = idx%nqz;
    const int qyi = (idx/nqz)%nqy;
    const real qxs=tex1Dfetch(tQx,qxi), qys=tex1Dfetch(tQy,qyi), qzs=tex1Dfetch(tQz,qzi);

    real Re = 0;
    real Im = 0;
    int densityidx = 0;
    for (int xi=0; xi < nx; xi++) {
        const real xQx = qxs*tex1Dfetch(tx,xi);
        for (int yi=0; yi < ny; yi++) {
#if 0 // 46 s (660 Gflops, 620 W) without loop unrolling
            const real yQy = qys*tex1Dfetch(ty,yi);
            for (int zi=0; zi < nz; zi++) {
                const real QdotR = xQx+yQy+qzs*tex1Dfetch(tz,zi);
                const real rho = tex1Dfetch(tdensity,densityidx);
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
                densityidx++;
            }
#else // 37 s (820 Gflops, 650 W) with loop unrolling
            const real xyQxy = xQx + qys*tex1Dfetch(ty,yi);
            int zi=0, tail = nz%4;
            for (zi=0; zi < nz-tail; zi+=4) {
                const real QdotR0 = xyQxy+qzs*tex1Dfetch(tz,zi);
                const real QdotR1 = xyQxy+qzs*tex1Dfetch(tz,zi+1);
                const real QdotR2 = xyQxy+qzs*tex1Dfetch(tz,zi+2);
                const real QdotR3 = xyQxy+qzs*tex1Dfetch(tz,zi+3);
                const real rho0 = tex1Dfetch(tdensity,densityidx);
                const real rho1 = tex1Dfetch(tdensity,densityidx+1);
                const real rho2 = tex1Dfetch(tdensity,densityidx+2);
                const real rho3 = tex1Dfetch(tdensity,densityidx+3);
                Re += rho0*cos(QdotR0) + rho1*cos(QdotR1) + rho2*cos(QdotR2) + rho3*cos(QdotR3);
                Im += rho0*sin(QdotR0) + rho1*sin(QdotR1) + rho2*sin(QdotR2) + rho3*sin(QdotR3);
                densityidx+=4;
            }
            for (;zi < nz; zi++) {
                const real QdotR = xyQxy+qzs*tex1Dfetch(tz,zi);
                const real rho = tex1Dfetch(tdensity,densityidx);
                Re += rho*cos(QdotR);
                Im += rho*sin(QdotR);
                densityidx++;
            }
#endif
        }
    }
    result[idx].x = Re;
    result[idx].y = Im;
}
