#include "head.h"

#define tpb 256

extern double *d_t;
extern double *d_it;
extern double *d_V;
extern double *d_dV2;
extern double *d_Vnew;
extern double *d_m;
extern double *d_h;
extern double *d_jj;
extern double *d_d;
extern double *d_f;
extern double *d_X;
extern double *d_cai;

__global__ void boundary(double *d_V){
	int k = blockDim.x * blockIdx.x + threadIdx.x;

	if(k<nx){

	d_V[(k+1)*(nx+2)] = d_V[(k+1)*(nx+2)+1];
        d_V[(k+1)*(nx+2)+(nx+1)] = d_V[(k+1)*(nx+2)+nx];
        d_V[k+1] = d_V[k+1+(nx+2)];
        d_V[(ny+1)*(nx+2)+k+1] = d_V[ny*(nx+2)+k+1];

	}
}

__device__ void comp_it(double *d_V, double *d_m, double *d_h, double *d_jj,
                                 double *d_d, double *d_f, double *d_cai,
                                         double *d_X, double *d_it, int k, int j, int i) {

	//int id = k+nx+2+1+2*j;

	d_it[k] = 0.0;

	double gna = 23.0;
        double ena = ((R*temp) / frdy)*log(nao / nai);

        double am = 0.32*(d_V[k+nx+2+1+2*j] + 47.13) / (1 - exp(-0.1*(d_V[k+nx+2+1+2*j] + 47.13)));
        double bm = 0.08*exp(-d_V[k+nx+2+1+2*j] / 11);
	double ah, bh, aj ,bj;
        if (d_V[k+nx+2+1+2*j] < -40.0) {
                ah = 0.135*exp((80 + d_V[k+nx+2+1+2*j]) / -6.8);
                bh = 3.56*exp(0.079*d_V[k+nx+2+1+2*j]) + 310000 * exp(0.35*d_V[k+nx+2+1+2*j]);
                aj = (-127140 * exp(0.2444*d_V[k+nx+2+1+2*j]) - 0.00003474*exp(-0.04391*d_V[k+nx+2+1+2*j]))*
                        ((d_V[k+nx+2+1+2*j] + 37.78)/(1 + exp(0.311*(d_V[k+nx+2+1+2*j] + 79.23))));
                bj = (0.1212*exp(-0.01052*d_V[k+nx+2+1+2*j])) / (1 + exp(-0.1378*(d_V[k+nx+2+1+2*j] + 40.14)));
        }
        else {
                ah = 0;
                bh = 1 / (0.13*(1 + exp((d_V[k+nx+2+1+2*j] + 10.66) / -11.1)));
                aj = 0;
                bj = (0.3*exp(-0.0000002535*d_V[k+nx+2+1+2*j])) / (1 + exp(-0.1*(d_V[k+nx+2+1+2*j] + 32)));
        }
        double mtau = 1 / (am + bm);
        double htau = 1 / (ah + bh);
	double jtau = 1 / (aj + bj);

        double mss = am*mtau;
        double hss = ah*htau;
        double jss = aj*jtau;

        d_m[k] = mss - (mss - d_m[k])*exp(-dt / mtau);
        d_h[k] = hss - (hss - d_h[k])*exp(-dt / htau);
        d_jj[k] = jss - (jss - d_jj[k])*exp(-dt / jtau);

        //d_it[k] += gna*d_m[k] * d_m[k] * d_m[k] * d_h[k] * d_jj[k] * (d_V[k+nx+2+1+2*j] - ena);

	__shared__ double esi[tpb];
	__shared__ double isi[tpb];
        esi[i] = 7.7 - 13.0287*log(d_cai[k]);

        double ad = 0.095*exp(-0.01*(d_V[k+nx+2+1+2*j] - 5)) / (1 + exp(-0.072*(d_V[k+nx+2+1+2*j] - 5)));
        double bd = 0.07*exp(-0.017*(d_V[k+nx+2+1+2*j] + 44)) / (1 + exp(0.05*(d_V[k+nx+2+1+2*j] + 44)));
        double af = 0.012*exp(-0.008*(d_V[k+nx+2+1+2*j] + 28)) / (1 + exp(0.15*(d_V[k+nx+2+1+2*j] + 28)));
        double bf = 0.0065*exp(-0.02*(d_V[k+nx+2+1+2*j] + 30)) / (1 + exp(-0.2*(d_V[k+nx+2+1+2*j] + 30)));

        double taud = 1 / (ad + bd);
        double tauf = 1 / (af + bf);

        double dss = ad*taud;
        double fss = af*tauf;

        d_d[k] = dss - (dss - d_d[k])*exp(-dt / taud);
        d_f[k] = fss - (fss - d_f[k])*exp(-dt / tauf);

        isi[i] = 0.09*d_d[k] * d_f[k] * (d_V[k+nx+2+1+2*j] - esi[i]);

        double dcai = -0.0001*isi[i] + 0.07*(0.0001 - d_cai[k]);

        d_cai[k] = d_cai[k] + dcai*dt;
	//d_it[k] = d_it[k] + isi[i];

        double gk = 0.282*sqrt(ko / 5.4);
        double ek = ((R*temp) / frdy)*log(ko / ki);
        //double prnak = 0.01833;
        //ek = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));

        double ax = 0.0005*exp(0.083*(d_V[k+nx+2+1+2*j] + 50)) / (1 + exp(0.057*(d_V[k+nx+2+1+2*j] + 50)));
        double bx = 0.0013*exp(-0.06*(d_V[k+nx+2+1+2*j] + 20)) / (1 + exp(-0.04*(d_V[k+nx+2+1+2*j] + 20)));

        double taux = 1 / (ax + bx);
        double xss = ax*taux;
        d_X[k] = xss - (xss - d_X[k])*exp(-dt / taux);

	double Xi;
        if (d_V[k+nx+2+1+2*j] > -100) {
                Xi = 2.837*(exp(0.04*(d_V[k+nx+2+1+2*j] + 77)) - 1)/
			((d_V[k+nx+2+1+2*j] + 77)*exp(0.04*(d_V[k+nx+2+1+2*j] + 35)));
        }
        else {
                Xi = 1;
        }
        //d_it[k] += gk*d_X[k] * Xi*(d_V[k+nx+2+1+2*j] - ek);

        double gk1 = 0.6047*(sqrt(ko / 5.4));
        double ek1 = ((R*temp) / frdy)*log(ko / ki);

        double ak1 = 1.02 / (1 + exp(0.2385*(d_V[k+nx+2+1+2*j] - ek1 - 59.215)));
        double bk1 = (0.49124*exp(0.08032*(d_V[k+nx+2+1+2*j] - ek1 + 5.476))+
			exp(0.06175*(d_V[k+nx+2+1+2*j] - ek1 - 594.31)))
                	/(1 + exp(-0.5143*(d_V[k+nx+2+1+2*j] - ek1 + 4.753)));
        double K1ss = ak1 / (ak1 + bk1);

        //d_it[k] += gk1*K1ss*(d_V[k+nx+2+1+2*j] - ek1);

        double gkp = 0.0183;
        double ekp = ((R*temp) / frdy)*log(ko / ki);

        double kp = 1 / (1 + exp((7.488 - d_V[k+nx+2+1+2*j]) / 5.98));

        //d_it[k] += gkp*kp*(d_V[k+nx+2+1+2*j] - ekp);

        //d_it[k] += 0.03921*(d_V[k+nx+2+1+2*j] + 59.87);
	d_it[k] = d_it[k] + gna*d_m[k] * d_m[k] * d_m[k] * d_h[k] * d_jj[k] * (d_V[k+nx+2+1+2*j] - ena)
			+ isi[i]
			+ gk*d_X[k] * Xi*(d_V[k+nx+2+1+2*j] - ek)
			+ gk1*K1ss*(d_V[k+nx+2+1+2*j] - ek1)
			+ gkp*kp*(d_V[k+nx+2+1+2*j] - ekp)
			+ 0.03921*(d_V[k+nx+2+1+2*j] + 59.87);
}

__global__ void comp(double *d_V, double *d_m, double *d_h, double *d_jj,
			 double *d_it, double *d_d, double *d_f, double *d_cai,
				double *d_X, double *d_dV2){
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	int i = threadIdx.x;

        if(k<nx*ny){

        int j = (int)(k/nx);

	comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, d_X, d_it, k, j, i);
	}
}

__global__ void comp_dV2(double *d_V ,double *d_dV2  ,double *d_it){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nx*ny){

	int j = (int)(k/nx);

	int id = k+(nx+2)+1+(2*j);

	d_dV2[k] = -d_it[k] + D*((d_V[id+1] + d_V[id-1] - 2*d_V[id])
                                 / (dx*dx) +(d_V[id+(nx+2)] + d_V[id-(nx+2)]-2*d_V[id])/(dy*dy));

	}
}

void gpu(){
	int bpg;
	//tpb = 256;
	bpg = (nx+tpb-1)/tpb;
	boundary<<<bpg, tpb>>>(d_V);
        bpg = (nx*ny+tpb-1)/tpb;
	comp<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_it, d_d, d_f, d_cai, d_X, d_dV2);
	comp_dV2<<<bpg, tpb>>>(d_V, d_dV2, d_it);

	cudaDeviceSynchronize();
}

__global__ void plane_waves(double *d_dV2){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<ny*5){
	int i, j;
	i = (int)(k/nx);
	j = k-i*nx;

	d_dV2[j*ny+i] = d_dV2[j*ny+i] + (-st);

	}
}

void stimu(){
	int bpg;
        //int tpb;

        //tpb = 256;
        bpg = (ny*5+tpb-1)/tpb;
	plane_waves<<<bpg, tpb>>>(d_dV2);
	cudaDeviceSynchronize();
}

__global__ void Euler(double *d_V, double *d_dV2, double *d_Vnew, double *d_t){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nx*ny){

	int j = (int)(k/nx);
	d_Vnew[k] = d_V[k+nx+2+1+2*j] + dt*d_dV2[k];
        d_V[k+nx+2+1+2*j] = d_Vnew[k];

	}

	if(k==0){

	d_t[0] = d_t[0] + dt;

	}
}

void Forward_Euler(){
	int bpg;
        //int tpb;

        //tpb = 256;
        bpg = (nx*ny+tpb-1)/tpb;
	Euler<<<bpg, tpb>>>(d_V, d_dV2, d_Vnew, d_t);
	cudaDeviceSynchronize();
}
