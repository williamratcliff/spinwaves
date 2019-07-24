__constant__ Real Qx[QXSIZE], Qy[QYSIZE], Qz[QZSIZE], x[XSIZE], y[YSIZE], z[ZSIZE];
extern "C" __global__ void
cudaBorn(const Real density[XSIZE][YSIZE][ZSIZE], int qxi, Cplx result[])
{
    const Cplx I(0,1);
    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= QYSIZE*QZSIZE) return;
    const int qzi = idx%QZSIZE;
    const int qyi = (idx/QZSIZE)%QYSIZE;

//printf("Entering kernel\n");
    Cplx sum(0,0);
    for (int xi=0; xi < XSIZE; xi++) {
        for (int yi=0; yi < YSIZE; yi++) {
#if 0
# pragma unroll 4
            for (int zi=0; zi < ZSIZE; zi++) {
                const Real QdotR = x[xi]*Qx[qxi]+y[yi]*Qy[qyi]+z[zi]*Qz[qzi];
                const Real rho = density[xi][yi][zi];
                //sum += Cplx(rho*cos(QdotR),rho*sin(QdotR));
                sum += rho*exp(I*QdotR);
            }
#else
            const Real xyQxy = x[xi]*Qx[qxi] + y[yi]*Qy[qyi];
            int zi;
            for (zi=0; zi < ZSIZE-(ZSIZE%4); zi+=4) {
                const Real QdotR0 = xyQxy+z[zi]*Qz[qzi];
                const Real QdotR1 = xyQxy+z[zi+1]*Qz[qzi];
                const Real QdotR2 = xyQxy+z[zi+2]*Qz[qzi];
                const Real QdotR3 = xyQxy+z[zi+3]*Qz[qzi];
                const Real rho0 = density[xi][yi][zi];
                const Real rho1 = density[xi][yi][zi+1];
                const Real rho2 = density[xi][yi][zi+2];
                const Real rho3 = density[xi][yi][zi+3];
                sum += Cplx(
                    rho0*cos(QdotR0) + rho1*cos(QdotR1) + rho2*cos(QdotR2) + rho3*cos(QdotR3),
                    rho0*sin(QdotR0) + rho1*sin(QdotR1) + rho2*sin(QdotR2) + rho3*sin(QdotR3));
            }
            for (; zi < ZSIZE; zi++) {
                const Real QdotR = xyQxy+z[zi]*Qz[qzi];
                const Real rho = density[xi][yi][zi];
                //sum += Cplx(rho*cos(QdotR),rho*sin(QdotR));
                sum += rho*exp(I*QdotR);
            }
#endif
        }
    }
    result[idx] = sum;
}
