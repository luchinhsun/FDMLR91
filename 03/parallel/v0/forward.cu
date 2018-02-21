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

extern double *d_m0;
extern double *d_h0;
extern double *d_jj0;
extern double *d_d0;
extern double *d_f0;
extern double *d_X0;

extern double *d_dVdt;
extern double *dcai;

extern double *d_belta;
extern double *d_ADIf;
extern double *d_y;

/*
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
*/
/*
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
*/
__device__ void comp_it(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *dcai, double *d_X, double *d_it, double *d_m0, double *d_h0, double *d_jj0, double *d_d0, double *d_f0, double *d_X0, int I, int i, int k, double *d_t) {
	//int id = k+nx+2+1+2*j;
	d_it[k] = 0.0;

	//comp_ina
	double gna = 23;
        double ena = ((R*temp) / frdy)*log(nao / nai);

        double am = 0.32*(d_V[k] + 47.13) / (1 - exp(-0.1*(d_V[k] + 47.13)));
        double bm = 0.08*exp(-d_V[k] / 11);
	double ah, bh, aj ,bj;
        if (d_V[k] < -40.0) {
                ah = 0.135*exp((80 + d_V[k]) / -6.8);
                bh = 3.56*exp(0.079*d_V[k]) + 310000 * exp(0.35*d_V[k]);
                aj = (-127140 * exp(0.2444*d_V[k]) - 0.00003474*exp(-0.04391*d_V[k]))*((d_V[k] + 37.78)/(1 + exp(0.311*(d_V[k] + 79.23))));
                bj = (0.1212*exp(-0.01052*d_V[k])) / (1 + exp(-0.1378*(d_V[k] + 40.14)));
        }
        else {
                ah = 0;
                bh = 1 / (0.13*(1 + exp((d_V[k] + 10.66) / -11.1)));
                aj = 0;
                bj = (0.3*exp(-0.0000002535*d_V[k])) / (1 + exp(-0.1*(d_V[k] + 32)));
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

        d_it[k] += gna*d_m0[k] * d_m0[k] * d_m0[k] * d_h0[k] * d_jj0[k] * (d_V[k] - ena);
	//comp_ical
	__shared__ double esi[tpb];
	__shared__ double isi[tpb];
        esi[I] = 7.7 - 13.0287*log(d_cai[k]);

        double ad = 50 * 0.095*exp(-0.01*(d_V[k] - 5)) / (1 + exp(-0.072*(d_V[k] - 5)));
        double bd = 50 * 0.07*exp(-0.017*(d_V[k] + 44)) / (1 + exp(0.05*(d_V[k] + 44)));
        double af = 50 * 0.012*exp(-0.008*(d_V[k] + 28)) / (1 + exp(0.15*(d_V[k] + 28)));
        double bf = 50 * 0.0065*exp(-0.02*(d_V[k] + 30)) / (1 + exp(-0.2*(d_V[k] + 30)));

        double taud = 1 / (ad + bd);
        double tauf = 1 / (af + bf);

        double dss = ad*taud;
        double fss = af*tauf;

        d_d0[k] = dss - (dss - d_d[k])*exp(-d_t[k] / taud);
        d_f0[k] = fss - (fss - d_f[k])*exp(-d_t[k] / tauf);

        isi[I] = 0.09*d_d0[k] * d_f0[k] * (d_V[k] - esi[I]);

        dcai[k] = -0.0001*isi[I] + 0.07*(0.0001 - d_cai[k]);

        //d_cai[k] = d_cai[k] + dcai*dt;
	d_it[k] = d_it[k] + isi[I];
	//comp_ik
        double gk = 0.282*sqrt(ko / 5.4);
        double ek = ((R*temp) / frdy)*log(ko / ki);
        //double prnak = 0.01833;
        //ek = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));

        double ax = 50 * 0.0005*exp(0.083*(d_V[k] + 50)) / (1 + exp(0.057*(d_V[k] + 50)));
        double bx = 50 * 0.0013*exp(-0.06*(d_V[k] + 20)) / (1 + exp(-0.04*(d_V[k] + 20)));

        double taux = 1 / (ax + bx);
        double xss = ax*taux;
        d_X0[k] = xss - (xss - d_X[k])*exp(-d_t[k] / taux);

	double Xi;
        if (d_V[k] > -100) {
                Xi = 2.837*(exp(0.04*(d_V[k] + 77)) - 1)/((d_V[k] + 77)*exp(0.04*(d_V[k] + 35)));
        }
        else {
                Xi = 1;
        }
        d_it[k] += gk*d_X0[k] * Xi*(d_V[k] - ek);
	//comp_ik1
        double gk1 = 0.6047*(sqrt(ko / 5.4));
        double ek1 = ((R*temp) / frdy)*log(ko / ki);

        double ak1 = 1.02 / (1 + exp(0.2385*(d_V[k] - ek1 - 59.215)));
        double bk1 = (0.49124*exp(0.08032*(d_V[k] - ek1 + 5.476))+exp(0.06175*(d_V[k] - ek1 - 594.31)))/(1 + exp(-0.5143*(d_V[k] - ek1 + 4.753)));
        double K1ss = ak1 / (ak1 + bk1);

        d_it[k] += gk1*K1ss*(d_V[k] - ek1);
	//comp_ikp
        double gkp = 0.0183;
        double ekp = ((R*temp) / frdy)*log(ko / ki);

        double kp = 1 / (1 + exp((7.488 - d_V[k]) / 5.98));

        d_it[k] += gkp*kp*(d_V[k] - ekp);
	//comp_ib
        d_it[k] += 0.03921*(d_V[k] + 59.87);

}

__global__ void comp_dVdt(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *dcai, double *d_X, double *d_it, double *d_m0, double *d_h0, double *d_jj0, double *d_d0, double *d_f0, double *d_X0, double *d_dVdt, double *d_t){

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


__device__ void gate(double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_X, double *d_m0, double *d_h0, double *d_jj0, double *d_d0, double *d_f0, double *d_X0, int k){
	d_m[k] = d_m0[k];
	d_h[k] = d_h0[k];
	d_jj[k] = d_jj0[k];
	d_d[k] = d_d0[k];
	d_f[k] = d_f0[k];
	d_X[k] = d_X0[k];
}

__global__ void comp_ODE_stim(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *dcai, double *d_X, double *d_it, double *d_m0, double *d_h0, double *d_jj0, double *d_d0, double *d_f0, double *d_X0, double *d_dVdt, double *d_t, int ncount){

	int k = threadIdx.x + blockIdx.x * blockDim.x;
        int I = threadIdx.x;

        if(k<nx*ny){

        int i = (int)(k/nx);
	int j = k - i*nx;
	int id = i*nx+j;
	int k1, k0, ttt;
	//int vid = (i+1)*(nx+2)+j+1;

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
		d_dVdt[id] = -d_it[id] + (-st)*(j>=0 && j<5 && ncount>=1 && ncount <= stimtime);
		d_V[id] = d_V[id] + d_t[id]*d_dVdt[id];
        }

	}
}

void ODE(int ncount){
	int bpg;
        bpg = (nx*ny+tpb-1)/tpb;
        comp_ODE_stim<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, d_dVdt, d_t, ncount);
	//bpg = ((nx-5)*ny+tpb-1)/tpb;
	//comp_ODE<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_m0, d_h0, d_jj0, d_d0, d_f0, d_X0, d_dVdt, d_t, 5);
}

__global__ void sweep_x(double *d_V, double *d_ADIf){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nx*ny){
	int i = (int)(k/nx);
	int j = k-nx*i;
	int id = i*nx+j;
	//d_ADIf[id] = 0.0;
	d_ADIf[id] = d_V[id]  + (dt_max*D / (dx*dx)/2)*((d_V[id] - 2 * d_V[id] + d_V[id+1])*(j==0)+(d_V[id-1] - 2 * d_V[id] + d_V[id])*(j==nx-1)+(d_V[id-1] - 2 * d_V[id] + d_V[id+1])*(j>0 && j<nx-1));
//	d_ADIf[k] = d_V[k]  + (dt_max*D / (dx*dx)/2)*((d_V[k] - 2 * d_V[k] + d_V[k+nx])*(k-nx*(int)(k/nx)==0)+(d_V[k-nx] - 2 * d_V[k] + d_V[k])*(k-nx*(int)(k/nx)==ny-1)+(d_V[k-nx] - 2 * d_V[k] + d_V[k+nx])*(k-nx*(int)(k/nx)>0 && k-nx*(int)(k/nx)<ny-1));
/*	if(j==0){
		d_ADIf[id] = d_V[id]  + (dt_max*D / (dx*dx)/2)*(d_V[id] - 2 * d_V[id] + d_V[id+1]);
	}else if(j==nx-1){
		d_ADIf[id] = d_V[id] + (dt_max*D / (dx*dx)/2)*(d_V[id-1] - 2 *d_V[id] + d_V[id]);
	}else{
		d_ADIf[id] = d_V[id] + (dt_max*D / (dx*dx)/2)*(d_V[id-1] - 2 * d_V[id] + d_V[id+1]);
	}*/

	}
/*
                for (j = 1; j < ny + 1; j++){
                        for (i = 1; i < nx + 1; i++){
                                if (j==1){
                                        f[i][j] = V[i][j]  + (eta/2)*(V[i][j] - 2 * V[i][j] + V[i][j + 1]);
                                }else if (j==ny){
                                        f[i][j] = V[i][j] + (eta/2)*(V[i][j - 1] - 2 * V[i][j] + V[i][j]);
                                }else{
                                        f[i][j] = V[i][j] + (eta/2)*(V[i][j - 1] - 2 * V[i][j] + V[i][j + 1]);
                                }
                        }
                }
*/
}

void x_direction(){

	int bpg;
        bpg = (nx*ny+tpb-1)/tpb;
	sweep_x<<<bpg, tpb>>>(d_V, d_ADIf);
}

__global__ void solve_x(double *d_V, double *d_ADIf, double *d_belta, double *d_y){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

        if(k<nx){
	double eta = dt_max*D / (dx*dx);
        double b = 1 + eta;
        double b_1 = 1 + eta / 2;//take care the boundary value
        double b_n = 1 + eta / 2;//take care the boundary value
        double c = -eta / 2;
        double a = -eta / 2;

	int i;

        d_belta[0+k] = c / b_1;//-dt_max*D / (2*dx*dx+dt_max*D);
        d_y[0+k] = d_ADIf[0+k] / b_1;//(1+dt_max*D / (2*dx*dx));
	for (i = 1; i < ny-1; i++){ //i = 2,3,...,n-1
        	d_belta[i*nx+k] = c/(b-a*d_belta[(i-1)*nx+k]);
        	d_y[i*nx+k] = (d_ADIf[i*nx+k] - a*d_y[(i-1)*nx+k]) / (b-a*d_belta[(i-1)*nx+k]);
        }
	d_y[(ny-1)*nx+k] = (d_ADIf[(ny-1)*nx+k] - a*d_y[(ny-1-1)*nx+k]) / (b_n - a*d_belta[(ny-1-1)*nx+k]);
        d_V[(ny-1)*nx+k] = d_y[(ny-1)*nx+k];
        for (i = ny-2; i >=0; i--){
        	d_V[i*nx+k] = d_y[i*nx+k] - d_belta[i*nx+k] * d_V[(i+1)*nx+k];
        }

	}

}

void x_solve(){

        int bpg;
        bpg = (nx+tpb-1)/tpb;
        solve_x<<<bpg, tpb>>>(d_V, d_ADIf, d_belta, d_y);
}

__global__ void sweep_y(double *d_V, double *d_ADIf){
        int k = threadIdx.x + blockIdx.x * blockDim.x;

        if(k<nx*ny){
        int i = (int)(k/nx);
        int j = k-nx*i;
        int id = i*nx+j;
        //d_ADIf[id] = 0.0;
        d_ADIf[id] = d_V[id]  + (dt_max*D / (dx*dx)/2)*((d_V[id] - 2 * d_V[id] + d_V[id+nx])*(i==0)+(d_V[id-nx] - 2 * d_V[id] + d_V[id])*(i==ny-1)+(d_V[id-nx] - 2 * d_V[id] + d_V[id+nx])*(i>0 && i<ny-1));
	}
}

void y_direction(){

        int bpg;
        bpg = (nx*ny+tpb-1)/tpb;
        sweep_y<<<bpg, tpb>>>(d_V, d_ADIf);
}

__global__ void solve_y(double *d_V, double *d_ADIf, double *d_belta, double *d_y){
        int k = threadIdx.x + blockIdx.x * blockDim.x;

        if(k<ny){
        double eta = dt_max*D / (dx*dx);
        double b = 1 + eta;
        double b_1 = 1 + eta / 2;//take care the boundary value
        double b_n = 1 + eta / 2;//take care the boundary value
        double c = -eta / 2;
        double a = -eta / 2;

        int j;

        d_belta[0+nx*k] = c / b_1;//-dt_max*D / (2*dx*dx+dt_max*D);
        d_y[0+nx*k] = d_ADIf[0+nx*k] / b_1;//(1+dt_max*D / (2*dx*dx));
        for (j = 1; j < nx-1; j++){ //i = 2,3,...,n-1
                d_belta[k*nx+j] = c/(b-a*d_belta[k*nx+j-1]);
                d_y[k*nx+j] = (d_ADIf[k*nx+j] - a*d_y[k*nx+j-1]) / (b-a*d_belta[k*nx+j-1]);
        }
        d_y[k*nx+nx-1] = (d_ADIf[k*nx+nx-1] - a*d_y[k*nx+nx-1-1]) / (b_n - a*d_belta[k*nx+nx-1-1]);
        d_V[k*nx+nx-1] = d_y[k*nx+nx-1];
        for (j = ny-2; j >=0; j--){
                d_V[k*nx+j] = d_y[k*nx+j] - d_belta[k*nx+j] * d_V[k*nx+j+1];
        }

        }

}

void y_solve(){

        int bpg;
        bpg = (ny+tpb-1)/tpb;
        solve_y<<<bpg, tpb>>>(d_V, d_ADIf, d_belta, d_y);
}
/*
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
*/
