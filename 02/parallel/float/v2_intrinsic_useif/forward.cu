#include "head.h"

#define tpb 256

extern float *d_t;
extern float *d_it;
extern float *d_V;
extern float *d_dV2;
extern float *d_Vnew;
extern float *d_m;
extern float *d_h;
extern float *d_jj;
extern float *d_d;
extern float *d_f;
extern float *d_X;
extern float *d_cai;

extern float *d_m0;
extern float *d_h0;
extern float *d_jj0;
extern float *d_d0;
extern float *d_f0;
extern float *d_X0;

extern float *d_dVdt;
extern float *dcai;

__global__ void boundary(float *d_V){
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


__global__ void comp_dV2(float *d_V ,float *d_dV2){
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

__device__ void comp_it(float *d_V, float *d_m, float *d_h, float *d_jj, float *d_d, float *d_f, float *d_cai, float *dcai, float *d_X, float *d_it, float *d_m0, float *d_h0, float *d_jj0, float *d_d0, float *d_f0, float *d_X0, int I, int i, int k, float *d_t) {
	//int id = k+nx+2+1+2*j;
	d_it[k] = 0.0;

	//comp_ina
	float gna = 23;
        float ena = ((R*temp) / frdy)*__logf(nao / nai);

        float am = 0.32*(d_V[k+nx+2+1+2*i] + 47.13) / (1 - __expf(-0.1*(d_V[k+nx+2+1+2*i] + 47.13)));
        float bm = 0.08*__expf(-d_V[k+nx+2+1+2*i] / 11);
	float ah, bh, aj ,bj;
        if (d_V[k+nx+2+1+2*i] < -40.0) {
                ah = 0.135*__expf((80 + d_V[k+nx+2+1+2*i]) / -6.8);
                bh = 3.56*__expf(0.079*d_V[k+nx+2+1+2*i]) + 310000 * __expf(0.35*d_V[k+nx+2+1+2*i]);
                aj = (-127140 * __expf(0.2444*d_V[k+nx+2+1+2*i]) - 0.00003474*__expf(-0.04391*d_V[k+nx+2+1+2*i]))*
                        ((d_V[k+nx+2+1+2*i] + 37.78)/(1 + __expf(0.311*(d_V[k+nx+2+1+2*i] + 79.23))));
                bj = (0.1212*__expf(-0.01052*d_V[k+nx+2+1+2*i])) / (1 + __expf(-0.1378*(d_V[k+nx+2+1+2*i] + 40.14)));
        }
        else {
                ah = 0;
                bh = 1 / (0.13*(1 + __expf((d_V[k+nx+2+1+2*i] + 10.66) / -11.1)));
                aj = 0;
                bj = (0.3*__expf(-0.0000002535*d_V[k+nx+2+1+2*i])) / (1 + __expf(-0.1*(d_V[k+nx+2+1+2*i] + 32)));
        }
        float mtau = 1 / (am + bm);
        float htau = 1 / (ah + bh);
	float jtau = 1 / (aj + bj);

        float mss = am*mtau;
        float hss = ah*htau;
        float jss = aj*jtau;

        d_m0[k] = mss - (mss - d_m[k])*__expf(-d_t[k] / mtau);
        d_h0[k] = hss - (hss - d_h[k])*__expf(-d_t[k] / htau);
        d_jj0[k] = jss - (jss - d_jj[k])*__expf(-d_t[k] / jtau);

        d_it[k] += gna*d_m0[k] * d_m0[k] * d_m0[k] * d_h0[k] * d_jj0[k] * (d_V[k+nx+2+1+2*i] - ena);
	//comp_ical
	__shared__ float esi[tpb];
	__shared__ float isi[tpb];
        esi[I] = 7.7 - 13.0287*__logf(d_cai[k]);

        float ad = 50 * 0.095*__expf(-0.01*(d_V[k+nx+2+1+2*i] - 5)) / (1 + __expf(-0.072*(d_V[k+nx+2+1+2*i] - 5)));
        float bd = 50 * 0.07*__expf(-0.017*(d_V[k+nx+2+1+2*i] + 44)) / (1 + __expf(0.05*(d_V[k+nx+2+1+2*i] + 44)));
        float af = 50 * 0.012*__expf(-0.008*(d_V[k+nx+2+1+2*i] + 28)) / (1 + __expf(0.15*(d_V[k+nx+2+1+2*i] + 28)));
        float bf = 50 * 0.0065*__expf(-0.02*(d_V[k+nx+2+1+2*i] + 30)) / (1 + __expf(-0.2*(d_V[k+nx+2+1+2*i] + 30)));

        float taud = 1 / (ad + bd);
        float tauf = 1 / (af + bf);

        float dss = ad*taud;
        float fss = af*tauf;

        d_d0[k] = dss - (dss - d_d[k])*__expf(-d_t[k] / taud);
        d_f0[k] = fss - (fss - d_f[k])*__expf(-d_t[k] / tauf);

        isi[I] = 0.09*d_d0[k] * d_f0[k] * (d_V[k+nx+2+1+2*i] - esi[I]);

        dcai[k] = -0.0001*isi[I] + 0.07*(0.0001 - d_cai[k]);

        //d_cai[k] = d_cai[k] + dcai*dt;
	d_it[k] = d_it[k] + isi[I];
	//comp_ik
        float gk = 0.282*sqrt(ko / 5.4);
        float ek = ((R*temp) / frdy)*__logf(ko / ki);
        //float prnak = 0.01833;
        //ek = ((R*temp) / frdy)*__logf((ko + prnak*nao) / (ki + prnak*nai));

        float ax = 50 * 0.0005*__expf(0.083*(d_V[k+nx+2+1+2*i] + 50)) / (1 + __expf(0.057*(d_V[k+nx+2+1+2*i] + 50)));
        float bx = 50 * 0.0013*__expf(-0.06*(d_V[k+nx+2+1+2*i] + 20)) / (1 + __expf(-0.04*(d_V[k+nx+2+1+2*i] + 20)));

        float taux = 1 / (ax + bx);
        float xss = ax*taux;
        d_X0[k] = xss - (xss - d_X[k])*__expf(-d_t[k] / taux);

	float Xi;
        if (d_V[k+nx+2+1+2*i] > -100) {
                Xi = 2.837*(__expf(0.04*(d_V[k+nx+2+1+2*i] + 77)) - 1)/((d_V[k+nx+2+1+2*i] + 77 + 1e-15)*__expf(0.04*(d_V[k+nx+2+1+2*i] + 35)));
        }
        else {
                Xi = 1;
        }
        d_it[k] += gk*d_X0[k] * Xi*(d_V[k+nx+2+1+2*i] - ek);
	//comp_ik1
        float gk1 = 0.6047*(sqrt(ko / 5.4));
        float ek1 = ((R*temp) / frdy)*__logf(ko / ki);

        float ak1 = 1.02 / (1 + __expf(0.2385*(d_V[k+nx+2+1+2*i] - ek1 - 59.215)));
        float bk1 = (0.49124*__expf(0.08032*(d_V[k+nx+2+1+2*i] - ek1 + 5.476))+__expf(0.06175*(d_V[k+nx+2+1+2*i] - ek1 - 594.31)))/(1 + __expf(-0.5143*(d_V[k+nx+2+1+2*i] - ek1 + 4.753)));
        float K1ss = ak1 / (ak1 + bk1);

        d_it[k] += gk1*K1ss*(d_V[k+nx+2+1+2*i] - ek1);
	//comp_ikp
        float gkp = 0.0183;
        float ekp = ((R*temp) / frdy)*__logf(ko / ki);

        float kp = 1 / (1 + __expf((7.488 - d_V[k+nx+2+1+2*i]) / 5.98));

        d_it[k] += gkp*kp*(d_V[k+nx+2+1+2*i] - ekp);
	//comp_ib
        d_it[k] += 0.03921*(d_V[k+nx+2+1+2*i] + 59.87);

}

__global__ void comp_dVdt(float *d_V, float *d_m, float *d_h, float *d_jj, float *d_d, float *d_f, float *d_cai, float *dcai, float *d_X, float *d_it, float *d_m0, float *d_h0, float *d_jj0, float *d_d0, float *d_f0, float *d_X0, float *d_dVdt, float *d_t){

	int k = threadIdx.x + blockIdx.x * blockDim.x;
        int I = threadIdx.x;

        if(k<nx*ny){

        int i = (int)(k/nx);
	d_t[k] = dt_max;
	comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, I, i, k, d_t);
	d_dVdt[k] = -d_it[k];
	}
}

void dVdt(){
	int bpg;

        bpg = (nx*ny+tpb-1)/tpb;
        comp_dVdt<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, d_dVdt, d_t);
}

__global__ void plane_waves(float *d_dVdt){
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


__device__ void gate(float *d_m, float *d_h, float *d_jj, float *d_d, float *d_f, float *d_X, float *d_m0, float *d_h0, float *d_jj0, float *d_d0, float *d_f0, float *d_X0, int k){
	d_m[k] = d_m0[k];
	d_h[k] = d_h0[k];
	d_jj[k] = d_jj0[k];
	d_d[k] = d_d0[k];
	d_f[k] = d_f0[k];
	d_X[k] = d_X0[k];
}

__global__ void comp_ODE_stim(float *d_V, float *d_m, float *d_h, float *d_jj, float *d_d, float *d_f, float *d_cai, float *dcai, float *d_X, float *d_it, float *d_m0, float *d_h0, float *d_jj0, float *d_d0, float *d_f0, float *d_X0, float *d_dVdt, float *d_t){

	int k = threadIdx.x + blockIdx.x * blockDim.x;
        int I = threadIdx.x;

        if(k<nx*ny){

        int i = (int)(k/nx);
	int j = k - i*nx;
	int id = i*nx+j;
	int k1, k0, ttt;
	int vid = (i+1)*(nx+2)+j+1;

	if(d_dVdt[id]>0){
		k0 = 5;
	}else{
		k0 = 1;
	}
	k1 = k0 + (int)(fabs(d_dVdt[id]) + 0.5);
	if (k1 >(int)(dt_max / dt_min)){
		k1 = (int)(dt_max / dt_min);
	}
	d_t[id] = dt_max / k1;
	for (ttt = 0; ttt < k1; ttt++){ //from t to t+dt_max, t=t+dt
		comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, I, i, id, d_t);
		gate(d_m, d_h, d_jj, d_d, d_f, d_X, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, id);
		d_cai[id] = d_cai[id] + dcai[id]*d_t[id];//renew Cai
		if(i>0 && i<5){
			d_dVdt[id] = -d_it[id] + (-st);
		}else{
			d_dVdt[id] = -d_it[id];
		}
		d_V[vid] = d_V[vid] + d_t[id]*d_dVdt[id];
        }

	}
}

__global__ void comp_ODE(float *d_V, float *d_m, float *d_h, float *d_jj, float *d_d, float *d_f, float *d_cai, float *dcai, float *d_X, float *d_it, float *d_m0, float *d_h0, float *d_jj0, float *d_d0, float *d_f0, float *d_X0, float *d_dVdt, float *d_t){

        int k = threadIdx.x + blockIdx.x * blockDim.x;
        int I = threadIdx.x;

        if(k<nx*ny){

        int i = (int)(k/nx);
        int j = k - i*nx;
        int id = i*nx+j;
        int k1, k0, ttt;
	int vid = (i+1)*(nx+2)+j+1;

        if(d_dVdt[id]>0){
                k0 = 5;
        }else{
                k0 = 1;
        }
        k1 = k0 + (int)(fabs(d_dVdt[id])+0.5);
        if (k1 >(int)(dt_max / dt_min)){
                k1 = (int)(dt_max / dt_min);
        }
        d_t[id] = dt_max / k1;
        for (ttt = 0; ttt < k1; ttt++){ //from t to t+dt_max, t=t+dt
                comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, I, i, id, d_t);
                gate(d_m, d_h, d_jj, d_d, d_f, d_X, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, id);
                d_cai[id] = d_cai[id] + dcai[id]*d_t[id];//renew Cai
                d_dVdt[id] = -d_it[id];
                d_V[vid] = d_V[vid] + d_t[id]*d_dVdt[id];
        }

        }
}

void ODE_stim(){
	int bpg;
        bpg = (nx*ny+tpb-1)/tpb;
        comp_ODE_stim<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, d_dVdt, d_t);
//	bpg = ((nx-5)*ny+tpb-1)/tpb;
//	comp_ODE<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, d_dVdt, d_t, 5);
}

void ODE(){
	int bpg;
        bpg = (nx*ny+tpb-1)/tpb;
        comp_ODE<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, d_dVdt, d_t);
}

__global__ void Euler(float *d_V, float *d_dV2, float *d_Vnew){
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
