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
/*
extern double *d_m0;
extern double *d_h0;
extern double *d_jj0;
extern double *d_d0;
extern double *d_f0;
extern double *d_X0;
*/
extern double *d_dVdt;
//extern double *dcai;
extern double *d_isi;
extern double *d_D1V;
extern double *d_D2V;

//extern int ncount;
__global__ void boundary(double *d_V){
	int k = blockDim.x * blockIdx.x + threadIdx.x;

	if(k<nx){

	d_V[(k+1)*(nx+2)] = d_V[(k+1)*(nx+2)+1];
        d_V[(k+1)*(nx+2)+(nx+1)] = d_V[(k+1)*(nx+2)+nx];
        d_V[k+1] = d_V[k+1+(nx+2)];
        d_V[(ny+1)*(nx+2)+k+1] = d_V[ny*(nx+2)+k+1];

	}
}

void bc(){
        int bpg;
        //tpb = 256;
        bpg = (nx+tpb-1)/tpb;
        boundary<<<bpg, tpb>>>(d_V);
        //cudaDeviceSynchronize();
}


__global__ void comp_dV2(double *d_V ,double *d_dV2){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nx*ny){

	int i = (int)(k/nx);
	int id = k+(nx+2)+1+(2*i);

	d_dV2[k] = D*((d_V[id+1] + d_V[id-1] - 2*d_V[id]) / (dx*dx) + (d_V[id+(nx+2)] + d_V[id-(nx+2)] - 2*d_V[id])/(dy*dy));

	}
}

void dV2(){
	int bpg;
	//tpb = 256;
        bpg = (nx*ny+tpb-1)/tpb;
	comp_dV2<<<bpg, tpb>>>(d_V, d_dV2);
	//cudaDeviceSynchronize();
}

__device__ void comp_it(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *d_isi, double *d_X, double *d_it, int I, int i, int k, double *d_t) {
	d_it[k] = 0.0;

	int id = k+nx+2+1+2*i;

	//comp_ina
	//double gna = 23;
        //double ena = ((R*temp) / frdy)*log(nao / nai);
	/*
        double am = 0.32*(d_V[k+nx+2+1+2*i] + 47.13) / (1 - exp(-0.1*(d_V[k+nx+2+1+2*i] + 47.13)));
        double bm = 0.08*exp(-d_V[k+nx+2+1+2*i] / 11);
	double ah, bh, aj ,bj;
        if (d_V[k+nx+2+1+2*i] < -40.0) {
                ah = 0.135*exp((80 + d_V[k+nx+2+1+2*i]) / -6.8);
                bh = 3.56*exp(0.079*d_V[k+nx+2+1+2*i]) + 310000 * exp(0.35*d_V[k+nx+2+1+2*i]);
                aj = (-127140 * exp(0.2444*d_V[k+nx+2+1+2*i]) - 0.00003474*exp(-0.04391*d_V[k+nx+2+1+2*i]))*
                        ((d_V[k+nx+2+1+2*i] + 37.78)/(1 + exp(0.311*(d_V[k+nx+2+1+2*i] + 79.23))));
                bj = (0.1212*exp(-0.01052*d_V[k+nx+2+1+2*i])) / (1 + exp(-0.1378*(d_V[k+nx+2+1+2*i] + 40.14)));
        }
        else {
                ah = 0;
                bh = 1 / (0.13*(1 + exp((d_V[k+nx+2+1+2*i] + 10.66) / -11.1)));
                aj = 0;
                bj = (0.3*exp(-0.0000002535*d_V[k+nx+2+1+2*i])) / (1 + exp(-0.1*(d_V[k+nx+2+1+2*i] + 32)));
        }
        double mtau = 1 / (am + bm);
        double htau = 1 / (ah + bh);
	double jtau = 1 / (aj + bj);

        double mss = am*mtau;
        double hss = ah*htau;
        double jss = aj*jtau;

        d_m0[k] = mss - (mss - d_m[k])*exp(-d_t[k] / mtau);
        d_h0[k] = hss - (hss - d_h[k])*exp(-d_t[k] / htau);
        d_jj0[k] = jss - (jss - d_jj[k])*exp(-d_t[k] / jtau);
	*/
        d_it[k] += gna*d_m[k] * d_m[k] * d_m[k] * d_h[k] * d_jj[k] * (d_V[id] - ena);
	//comp_ical
	__shared__ double esi[tpb];
	//__shared__ double isi[tpb];
        esi[I] = 7.7 - 13.0287*log(d_cai[k]);
	/*
        double ad = 50 * 0.095*exp(-0.01*(d_V[k+nx+2+1+2*i] - 5)) / (1 + exp(-0.072*(d_V[k+nx+2+1+2*i] - 5)));
        double bd = 50 * 0.07*exp(-0.017*(d_V[k+nx+2+1+2*i] + 44)) / (1 + exp(0.05*(d_V[k+nx+2+1+2*i] + 44)));
        double af = 50 * 0.012*exp(-0.008*(d_V[k+nx+2+1+2*i] + 28)) / (1 + exp(0.15*(d_V[k+nx+2+1+2*i] + 28)));
        double bf = 50 * 0.0065*exp(-0.02*(d_V[k+nx+2+1+2*i] + 30)) / (1 + exp(-0.2*(d_V[k+nx+2+1+2*i] + 30)));

        double taud = 1 / (ad + bd);
        double tauf = 1 / (af + bf);

        double dss = ad*taud;
        double fss = af*tauf;

        d_d0[k] = dss - (dss - d_d[k])*exp(-d_t[k] / taud);
        d_f0[k] = fss - (fss - d_f[k])*exp(-d_t[k] / tauf);
	*/
        d_isi[k] = 0.09*d_d[k] * d_f[k] * (d_V[id] - esi[I]);

        //dcai[k] = -0.0001*isi[I] + 0.07*(0.0001 - d_cai[k]);

        //d_cai[k] = d_cai[k] + dcai*dt;
	d_it[k] = d_it[k] + d_isi[k];
	//comp_ik
	/*
        double gk = 0.282*sqrt(ko / 5.4);
        double ek = ((R*temp) / frdy)*log(ko / ki);
        //double prnak = 0.01833;
        //ek = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));

        double ax = 50 * 0.0005*exp(0.083*(d_V[k+nx+2+1+2*i] + 50)) / (1 + exp(0.057*(d_V[k+nx+2+1+2*i] + 50)));
        double bx = 50 * 0.0013*exp(-0.06*(d_V[k+nx+2+1+2*i] + 20)) / (1 + exp(-0.04*(d_V[k+nx+2+1+2*i] + 20)));

        double taux = 1 / (ax + bx);
        double xss = ax*taux;
        d_X0[k] = xss - (xss - d_X[k])*exp(-d_t[k] / taux);
	*/
	double Xi;
        if (d_V[id] > -100) {
                Xi = 2.837*(exp(0.04*(d_V[id] + 77)) - 1)/((d_V[id] + 77)*exp(0.04*(d_V[id] + 35)));
        }
        else {
                Xi = 1;
        }
        d_it[k] += gk*d_X[k] * Xi*(d_V[id] - ek);
	//comp_ik1
        //double gk1 = 0.6047*(sqrt(ko / 5.4));
        //double ek1 = ((R*temp) / frdy)*log(ko / ki);

        double ak1 = 1.02 / (1 + exp(0.2385*(d_V[id] - ek1 - 59.215)));
        double bk1 = (0.49124*exp(0.08032*(d_V[id] - ek1 + 5.476))+exp(0.06175*(d_V[id] - ek1 - 594.31)))/(1 + exp(-0.5143*(d_V[id] - ek1 + 4.753)));
        double K1ss = ak1 / (ak1 + bk1);

        d_it[k] += gk1*K1ss*(d_V[id] - ek1);
	//comp_ikp
        //double gkp = 0.0183;
        //double ekp = ((R*temp) / frdy)*log(ko / ki);

        double kp = 1 / (1 + exp((7.488 - d_V[id]) / 5.98));

        d_it[k] += gkp*kp*(d_V[id] - ekp);
	//comp_ib
        d_it[k] += 0.03921*(d_V[id] + 59.87);

}

__global__ void comp_dVdt(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *d_isi, double *d_X, double *d_it, double *d_dVdt, double *d_t){
	int k = threadIdx.x + blockIdx.x * blockDim.x;
        int I = threadIdx.x;

        if(k<nx*ny){

        int i = (int)(k/nx);
	d_t[k] = dt_max;
	comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, d_isi, d_X, d_it, I, i, k, d_t);
	d_dVdt[k] = -d_it[k];
	}
}

void dVdt(){
	int bpg;

        bpg = (nx*ny+tpb-1)/tpb;
        comp_dVdt<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, d_isi, d_X, d_it, d_dVdt, d_t);
}

__global__ void plane_waves(double *d_dVdt){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<ny*5){
	int i, j, id;
	i = (int)(k/5);
	j = k-i*5;
	id = i*nx+j;

	d_dVdt[id] = d_dVdt[id] + (-st);

	}
}

void stimu(){
	int bpg;
        //int tpb;

        //tpb = 256;
        bpg = (ny*5+tpb-1)/tpb;
	plane_waves<<<bpg, tpb>>>(d_dVdt);
	//cudaDeviceSynchronize();
}

/* Calculating the initial Ions Current*/
__device__ void comp_fluxes(double *d_cai, double *d_isi, int k, double *d_t){
        d_cai[k] = d_cai[k] + (-0.0001*d_isi[k] + 0.07*(0.0001 - d_cai[k]))*d_t[k];
}

/* update the gate value*/
__device__ void Rush_Larsen(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_X, int i, int k, double *d_t){
	//Fast sodium current
        //gate variables can not be shared, should be local due to data racing !!!!!!!!
	int id = k+nx+2+1+2*i;
	//double bm = 0.08*exp(-d_V[k+nx+2+1+2*i] / 11);
	double am = 0.32*(d_V[id] + 47.13) / (1 - exp(-0.1*(d_V[id] + 47.13)));
        double bm = 0.08*exp(-d_V[id] / 11);
        double ah, bh, aj, bj;
        if (d_V[id] < -40.0) {
                ah = 0.135*exp((80 + d_V[id]) / -6.8);
                bh = 3.56*exp(0.079*d_V[id]) + 310000.0 * exp(0.35*d_V[id]);
                aj = (-127140 * exp(0.2444*d_V[id]) - 0.00003474*exp(-0.04391*d_V[id]))*((d_V[id] + 37.78) / (1.0 + exp(0.311*(d_V[id] + 79.23))));
                bj = (0.1212*exp(-0.01052*d_V[id])) / (1.0 + exp(-0.1378*(d_V[id] + 40.14)));
        }
        else {
                ah = 0.0;
                bh = 1.0 / (0.13*(1.0 + exp((d_V[id] + 10.66) / -11.1)));
		aj = 0.0;
                bj = (0.3*exp(-0.0000002535*d_V[id])) / (1.0 + exp(-0.1*(d_V[id] + 32.0)));
        }
        double mtau = 1.0 / (am + bm);
        double htau = 1.0 / (ah + bh);
        double jtau = 1.0 / (aj + bj);

        double mss = am*mtau;
        double hss = ah*htau;
        double jss = aj*jtau;
        d_m[k] = mss - (mss - d_m[k])*exp(-d_t[k] / mtau);
	d_h[k] = hss - (hss - d_h[k])*exp(-d_t[k] / htau);
        d_jj[k] = jss - (jss - d_jj[k])*exp(-d_t[k] / jtau);


        //Slow inward current
        double ad = 50*0.095*exp(-0.01*(d_V[id] - 5)) / (1.0 + exp(-0.072*(d_V[id] - 5)));
        double bd = 50*0.07*exp(-0.017*(d_V[id] + 44)) / (1.0 + exp(0.05*(d_V[id] + 44)));
        double af = 50*0.012*exp(-0.008*(d_V[id] + 28)) / (1.0 + exp(0.15*(d_V[id] + 28)));
        double bf = 50*0.0065*exp(-0.02*(d_V[id] + 30)) / (1.0 + exp(-0.2*(d_V[id] + 30)));

        double taud = 1.0 / (ad + bd);
	double tauf = 1.0 / (af + bf);

        double dss = ad*taud;
        double fss = af*tauf;
        d_d[k] = dss - (dss - d_d[k])*exp(-d_t[k] / taud);
        d_f[k] = fss - (fss - d_f[k])*exp(-d_t[k] / tauf);

        //Time-dependent potassium current
        double ax = 50*0.0005*exp(0.083*(d_V[id] + 50)) / (1 + exp(0.057*(d_V[id] + 50)));
        double bx = 50*0.0013*exp(-0.06*(d_V[id] + 20)) / (1 + exp(-0.04*(d_V[id] + 20)));
        double taux = 1 / (ax + bx);
	double xss = ax*taux;
        d_X[k] = xss - (xss - d_X[k])*exp(-d_t[k] / taux);
}

__global__ void firsttime(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_X, double *d_cai, double *d_isi, double *d_dVdt, double *d_D1V, double *d_t){
	// in order to get D1V[i][j], for computing D2V[i][j] in CCL(i, j, dt_max);
	/* The first time step*/
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	int i = (int)(k/nx);
	int id = k+nx+2+1+2*i;

	if(k<nx*ny){
		d_D1V[k] = d_dVdt[k];
		comp_fluxes(d_cai, d_isi, k, d_t);
        	Rush_Larsen(d_V, d_m, d_h, d_jj, d_d, d_f, d_X, i, k, d_t);
        	d_V[id] = d_V[id] + dt_max * d_dVdt[k];
	}
}

void First(){
        int bpg;

        bpg = (nx*ny+tpb-1)/tpb;
        firsttime<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_X, d_cai, d_isi, d_dVdt, d_D1V, d_t);
        //cudaDeviceSynchronize();
}

__device__ void CCL(double *d_dVdt, double *d_D1V, double *d_D2V, int k, double *d_t){
/*
	double dt_range;
	dt_range = d_t[k]*2*(dt_univ > d_t[k]*2) + dt_univ*(dt_univ <= d_t[k]*2);

	d_D2V[k] = (d_dVdt[k] - d_D1V[k]) / d_t[k];
        double DiscriminantP = 0.0, DiscriminantN = 0.0, dtz = 0.0;
	DiscriminantP = d_dVdt[k] * d_dVdt[k] + 2 * d_D2V[k] * Voffset*(d_dVdt[k] >= 0);
        DiscriminantN = d_dVdt[k] * d_dVdt[k] - 2 * d_D2V[k] * Voffset*(d_dVdt[k] < 0);
        dtz = -d_dVdt[k] / d_D2V[k];
        d_t[k] = (-d_dVdt[k] + sqrt(DiscriminantP)) / d_D2V[k]*(d_dVdt[k] >= 0)*(d_D2V[k]>0)+(-d_dVdt[k] + sqrt(DiscriminantP)) / d_D2V[k]*(d_dVdt[k] >= 0)*(d_D2V[k]<0)*(DiscriminantP >= 0)+dtz*(d_dVdt[k] >= 0)*(d_D2V[k]<0)*(DiscriminantP < 0)+(-d_dVdt[k]+sqrt(DiscriminantN)) / d_D2V[k]*(d_dVdt[k] < 0)*(d_D2V[k]>0)*(DiscriminantN >= 0)+dtz*(d_dVdt[k] < 0)*(d_D2V[k]>0)*(DiscriminantN < 0)+(-d_dVdt[k] + sqrt(DiscriminantN)) / d_D2V[k]*(d_dVdt[k] < 0)*(d_D2V[k]<0);

        d_t[k] = d_t[k]*(d_t[k]<=dt_range && d_t[k]>=dt_min)+dt_range*(d_t[k]>dt_range)+dt_min*(d_t[k]<dt_min);
*/
	double dt_range;
        if (dt_univ > d_t[k] * 2){
                dt_range = d_t[k] * 2;
        }
        else{
                dt_range = dt_univ;
        }

        d_D2V[k] = (d_dVdt[k] - d_D1V[k]) / d_t[k];
        double DiscriminantP = 0, DiscriminantN = 0, dtz = 0;
        if (d_dVdt[k] >= 0){
                DiscriminantP = d_dVdt[k] * d_dVdt[k] + 2 * d_D2V[k] * Voffset;
                if (d_D2V[k]>0){
                        d_t[k] = (-d_dVdt[k] + sqrt(DiscriminantP)) / d_D2V[k];
                }
                else if (d_D2V[k]<0){
                        dtz = -d_dVdt[k] / d_D2V[k];
                        if (DiscriminantP >= 0){
                                d_t[k] = (-d_dVdt[k] + sqrt(DiscriminantP)) / d_D2V[k];
                        }
                        else{
                                d_t[k] = dtz;
                        }
                }
        }
	else{
                DiscriminantN = d_dVdt[k] * d_dVdt[k] - 2 * d_D2V[k] * Voffset;
                if (d_D2V[k]>0){
                        dtz = -d_dVdt[k] / d_D2V[k];
                        if (DiscriminantN >= 0){
                                d_t[k] = (-d_dVdt[k] - sqrt(DiscriminantN)) / d_D2V[k];
                        }
                        else{
                                d_t[k] = dtz;
                        }
                }
	else if (d_D2V[k]<0){
                        d_t[k] = (-d_dVdt[k] - sqrt(DiscriminantN)) / d_D2V[k];
                }
        }

        if (d_t[k]>dt_range){
                d_t[k] = dt_range;
        }
        if (d_t[k]<dt_min){
                d_t[k] = dt_min;
        }

}

__device__ void CCL_dtmax(double *d_dVdt, double *d_D1V, double *d_D2V, int k, double *d_t){
	double dt_range;
        if (dt_univ > dt_max * 2){
                dt_range = dt_max * 2;
        }
        else{
                dt_range = dt_univ;
        }

        d_D2V[k] = (d_dVdt[k] - d_D1V[k]) / dt_max;
        double DiscriminantP = 0, DiscriminantN = 0, dtz = 0;
        if (d_dVdt[k] >= 0){
                DiscriminantP = d_dVdt[k] * d_dVdt[k] + 2 * d_D2V[k] * Voffset;
                if (d_D2V[k]>0){
                        d_t[k] = (-d_dVdt[k] + sqrt(DiscriminantP)) / d_D2V[k];
                }
                else if (d_D2V[k]<0){
                        dtz = -d_dVdt[k] / d_D2V[k];
                        if (DiscriminantP >= 0){
                                d_t[k] = (-d_dVdt[k] + sqrt(DiscriminantP)) / d_D2V[k];
                        }
                        else{
                                d_t[k] = dtz;
                        }
                }
        }
	else{
                DiscriminantN = d_dVdt[k] * d_dVdt[k] - 2 * d_D2V[k] * Voffset;
                if (d_D2V[k]>0){
                        dtz = -d_dVdt[k] / d_D2V[k];
                        if (DiscriminantN >= 0){
                                d_t[k] = (-d_dVdt[k] - sqrt(DiscriminantN)) / d_D2V[k];
                        }
                        else{
                                d_t[k] = dtz;
                        }
                }
	else if (d_D2V[k]<0){
                        d_t[k] = (-d_dVdt[k] - sqrt(DiscriminantN)) / d_D2V[k];
                }
        }

        if (d_t[k]>dt_range){
                d_t[k] = dt_range;
        }
        if (d_t[k]<dt_min){
                d_t[k] = dt_min;
        }

}

__global__ void ODE_CCL(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *d_isi, double *d_X, double *d_it, double *d_dVdt, double *d_D1V, double *d_D2V, double *d_t, int ncount){
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	int i = (int)(k/nx);
	int j = k - i*nx;
	int id = k+nx+2+1+2*i;
	int I = threadIdx.x;

	double dt_sum;
	if(k<nx*ny){
		//***** adjust or correct time step---CCL method **/
                dt_sum = 0.0;
                CCL_dtmax(d_dVdt, d_D1V, d_D2V, k, d_t);
                dt_sum = dt_sum + d_t[k];
                while(dt_sum<dt_max){
                	d_D1V[k] = d_dVdt[k];
                	comp_fluxes(d_cai, d_isi, k, d_t);
                	Rush_Larsen(d_V, d_m, d_h, d_jj, d_d, d_f, d_X, i, k, d_t);
                	d_V[id] = d_V[id] + d_t[k] * d_dVdt[k] + d_t[k] * d_t[k] * d_D2V[k] / 2;
                	comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, d_isi, d_X, d_it, I, i, k, d_t);
                	d_dVdt[k] = -d_it[k] + (-st)*(ncount >= 1 && ncount <= stimtime && j >= 0 && j <= 4);
                	CCL(d_dVdt, d_D1V, d_D2V, k, d_t);
			dt_sum = dt_sum + d_t[k];
                }
                d_t[k] = dt_max - (dt_sum - d_t[k]);// here is a new dt  !!!
                d_D1V[k] = d_dVdt[k];
                comp_fluxes(d_cai, d_isi, k, d_t);
                Rush_Larsen(d_V, d_m, d_h, d_jj, d_d, d_f, d_X, i, k, d_t);
                d_V[id] = d_V[id] + d_t[k] * d_dVdt[k] + d_t[k] * d_t[k] * d_D2V[k] / 2;
	}
}

void ODE(int ncount){
        int bpg;

        bpg = (nx*ny+tpb-1)/tpb;
        ODE_CCL<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, d_isi, d_X, d_it, d_dVdt, d_D1V, d_D2V, d_t, ncount);
        //cudaDeviceSynchronize();
}

__global__ void Euler(double *d_V, double *d_dV2, double *d_Vnew){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nx*ny){

	int i = (int)(k/nx);
	d_Vnew[k] = d_V[k+nx+2+1+2*i] + dt_max/2 *d_dV2[k];
        d_V[k+nx+2+1+2*i] = d_Vnew[k];

	}
}

void Forward_Euler(){
	int bpg;
        //int tpb;

        //tpb = 256;
        bpg = (nx*ny+tpb-1)/tpb;
	Euler<<<bpg, tpb>>>(d_V, d_dV2, d_Vnew);
	//cudaDeviceSynchronize();
}
