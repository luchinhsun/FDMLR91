#include "head.h"

extern double *h_t;
extern double *h_V;
extern double *h_Vnew;
extern double *h_m;
extern double *h_h;
extern double *h_jj;
extern double *h_d;
extern double *h_f;
extern double *h_X;
extern double *h_cai;

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

omp_set_num_threads(2);

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

	int nstep = 500; // snapshot interval to save data files 500*0.02=10 ms
	int index = 0;// filename index from 1-5

	char filename[100];
	FILE *a;
	int fileflag;

	struct timeb start, end;
        int diff;
	ftime(&start);

	Send_to_Device();
	ncount = 0;
#pragma omp parallel private(i, j)
{

int tid = omp_get_thread_num();
	//for (ncount = 0; ncount <= 50000; ncount++){//30000 steps, 600ms
	while(ncount>=0 && ncount <= 50000){
	if(tid == 0){
		gpu();

		if (ncount >= 0 && ncount <= 100) {
                        stimu();
                }

		Forward_Euler();

		if (ncount == 20000){
			Send_to_Host();
                        for (i = 0; i < nx/2-1; i++){
                                for (j = 0; j < ny-1; j++){
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
			Send_to_Device();
                }
	ncount = ncount + 1;
	//if(ncount>50000)        break;
	}

	if(ncount%nstep == 0){//get data at the 9000th step
		#pragma omp barrier
		if(tid == 0){
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
		}
	}
//if(ncount>50000)	break;
	//printf("id = %d, ncount = %d\n", tid, ncount);
	}

}

	Send_V();
	ftime(&end);
        diff = (int)(1000.0*(end.time-start.time)+(end.millitm-start.millitm));
        printf("\nTime = %d ms\n", diff);
	Save_Result();
	free();
	return 0;
	//double time_used = (double)(end - start) / CLK_TCK;
	//fprintf(fevaluation, "%g", time_used);
}
