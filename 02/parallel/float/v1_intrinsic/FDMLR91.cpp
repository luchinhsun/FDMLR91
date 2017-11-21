#include "head.h"

extern float *h_t;
extern float *h_V;
extern float *h_Vnew;
extern float *h_m;
extern float *h_h;
extern float *h_jj;
extern float *h_d;
extern float *h_f;
extern float *h_X;
extern float *h_cai;

int main(){
	// Time Loop Conditions
	//h_t[0] = 0.0; // Time (ms)
	//	steps = (bcl*beats)/udt; // Number of ms
	//st = -80.0; // Stimulus (mA)

	// Beginning Ion Concentrations
	//nai = 18; // Initial Intracellular Na (mM)
	//nao = 140; // Initial Extracellular Na (mM)
	//ki = 145; // Initial Intracellular K (mM)
	//ko = 5.4; // Initial Extracellular K (mM)
	//cai = 0.0002; // Initial Intracellular Ca (mM)
	//cao = 1.8; // Initial Extracellular Ca (mM)

	int ncount, i, j;

	Allocate();
	h_t[0] = 0.0;
	h_V[0] = 0.0;h_V[nx+1] = 0.0;
	h_V[(nx+2)*(ny+1)] = 0.0;h_V[(nx+2)*(ny+1)+nx+1] = 0.0;
	for (i = 1; i < nx + 1; i++){
                h_V[i] = 0.0;
                h_V[i+(nx+2)*(ny+1)] = 0.0;
		h_V[0+(nx+2)*i] = 0.0;
                h_V[nx+1+i*(nx+2)] = 0.0;
        }

	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++){
			h_V[(i+1)*(nx+2)+j+1] = -88.654973; // Initial Voltage (mv)
			h_m[i*nx+j] = 0.000838;
			h_h[i*nx+j] = 0.993336;
			h_jj[i*nx+j] = 0.995484;
			h_d[i*nx+j] = 0.000003;
			h_f[i*nx+j] = 0.999745;
			h_X[i*nx+j] = 0.000129;
			h_cai[i*nx+j] = 0.0002; // Initial Intracellular Ca (mM)
		}
	}

	int nstep = 4/dt_max; // snapshot interval to save data files 500*0.02=10 ms
	int index = 0;// filename index from 1-5

	char filename[100];
	FILE *a;
	int fileflag;

	struct timeb start, end;
        int diff;
	ftime(&start);

	Send_to_Device();

	for (ncount = 0; ncount <= 160 / dt_max - 1; ncount++){//simulation time is 160ms
		bc();

		//performance();
/*		if (ncount%nstep == 0){//get data at the 9000th step
			Send_V();
			fileflag = 0;
			for (i = 0; i < nx; i++){
        			for (j = 0; j < ny; j++){
					if (fileflag == 0){
               	         			sprintf(filename, "a%d.txt", index);
                                		a = fopen(filename, "w");
                                		fileflag = 1;
                                		index++;
                        		}
                        		fprintf(a, "%g\t", h_Vnew[i*nx+j]);
                        		if (j == nx-1){
                                	        fprintf(a, "\n");
                        		}
				}
			}
			if (fileflag == 1){
        			fclose(a);
        		}
		}*/

		//-*********** step 1 *******
		dV2();
		Forward_Euler();

		//-*********** step 2 *******
		dVdt();
		if (ncount >= 1 && ncount <= stimtime) {
                        stimu();
                }
		if (ncount >= 1 && ncount <= stimtime) {
			ODE_stim();
		}else{
			ODE();
		}

		//-*********** step 3 *******
		dV2();
		Forward_Euler();

	}
	//Send_to_Host();
	Send_V();
	ftime(&end);
        diff = (int)(1000.0*(end.time-start.time)+(end.millitm-start.millitm));
        printf("\nTime = %d ms\n", diff);
	Save_Result();
	free();
	return 0;
	//float time_used = (float)(end - start) / CLK_TCK;
	//fprintf(fevaluation, "%g", time_used);
}
