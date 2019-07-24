
__global__ void

cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz,
		const real density[],
		const real x[], const real y[], const real z[],
		const real Qx[], const real Qy[], const real Qz[],
		int qxi,
		cplx result[])
{


    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    const int qzi = idx%nqz;
    const int qyi = (idx/nqz)%nqy;

    real precx,presx;
    real precy,presy;
    real precz,presz;

    real exx,exy,exz;

    exx = Qx[qxi]*(x[1]-x[0]);
    exy = Qy[qyi]*(y[1]-y[0]);
    exz = Qz[qzi]*(z[1]-z[0]);

    sincos(exx,&presx,&precx);
    sincos(exy,&presy,&precy);
    sincos(exz,&presz,&precz);

    real x_compRe = presx/Qx[qxi];
    real x_compIm = (1/Qx[qxi])*(1-precx);

    real y_compRe = presy/Qy[qyi];
	real y_compIm = (1/Qy[qyi])*(1-precy);

	real z_compRe = presz/Qz[qzi];
	real z_compIm = (1/Qz[qzi])*(1-precz);

	real compxyRe = ((x_compRe*y_compRe) - (x_compIm*y_compIm));
	real compxyIm = ((x_compIm*y_compRe) + (x_compRe*y_compIm));

	real totcompRe = ((compxyRe*z_compRe) - (compxyIm*z_compIm));
	real totcompIm = ((compxyIm*z_compRe) + (compxyRe*z_compIm));

	int xi, yi, zi;
    int densityidx = 0;

    real Re = 0;
    real Im = 0;


    if (qyi >= nqy || qzi >= nqz) return;

    for (xi=0; xi < nx; xi++) {
        for (yi=0; yi < ny; yi++) {
            for (zi=0; zi < nz; zi++) {

                const real QdotR = Qx[qxi]*x[xi] +Qy[qyi]*y[yi] + Qz[qzi]*z[zi];

                real cx,sx;
                real singRe,singIm;

                sincos(QdotR,&sx,&cx);

                singRe = density[densityidx]*cx;
                singIm = density[densityidx]*sx;

                Re += ((singRe*totcompRe)-(singIm*totcompIm));
                Im += ((singIm*totcompRe)+(singRe*totcompIm));
                densityidx++;
            }
        }
    }

    result[idx].x = Re;
    result[idx].y = Im;
}

