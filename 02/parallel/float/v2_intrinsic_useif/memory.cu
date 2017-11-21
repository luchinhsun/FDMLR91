#include "head.h"

float *h_t;
float *d_t;
float *h_V;
float *d_V;
float *d_dV2;
float *h_Vnew;
float *d_Vnew;
float *d_it;

float *h_m;
float *d_m;
float *h_h;
float *d_h;
float *h_jj;
float *d_jj;
float *h_d;
float *d_d;
float *h_f;
float *d_f;
float *h_X;
float *d_X;
float *h_cai;
float *d_cai;

float *h_it;

float *d_m0;
float *d_h0;
float *d_jj0;
float *d_d0;
float *d_f0;
float *d_X0;

float *d_dVdt;
float *dcai;

void Allocate(){
	cudaError_t Error;
	size_t size = nx*ny*sizeof(float);

	h_t = (float*)malloc(size);
	Error = cudaMalloc((void**)&d_t, size);
	printf("CUDA error = %s\n",cudaGetErrorString(Error));

	h_V = (float*)malloc((nx+2)*(ny+2)*sizeof(float));
	cudaMalloc((void**)&d_V, (nx+2)*(ny+2)*sizeof(float));
	cudaMalloc((void**)&d_dV2, size);
	h_Vnew = (float*)malloc(size);
	cudaMalloc((void**)&d_Vnew, size);

	cudaMalloc((void**)&d_it, size);

	h_m = (float*)malloc(size);
	cudaMalloc((void**)&d_m, size);
	h_h = (float*)malloc(size);
        cudaMalloc((void**)&d_h, size);
	h_jj = (float*)malloc(size);
        cudaMalloc((void**)&d_jj, size);
	h_d = (float*)malloc(size);
        cudaMalloc((void**)&d_d, size);
	h_f = (float*)malloc(size);
        cudaMalloc((void**)&d_f, size);
	h_X = (float*)malloc(size);
        cudaMalloc((void**)&d_X, size);
	h_cai = (float*)malloc(size);
        cudaMalloc((void**)&d_cai, size);

	h_it = (float*)malloc(size);

	cudaMalloc((void**)&d_m0, size);
        cudaMalloc((void**)&d_h0, size);
        cudaMalloc((void**)&d_jj0, size);
        cudaMalloc((void**)&d_d0, size);
        cudaMalloc((void**)&d_f0, size);
        cudaMalloc((void**)&d_X0, size);

	cudaMalloc((void**)&d_dVdt, size);
	cudaMalloc((void**)&dcai, size);
}

void free(){

	free(h_t);free(h_V);free(h_m);free(h_h);
	free(h_jj);free(h_d);free(h_f);free(h_X);free(h_cai);
	free(h_Vnew);
	free(h_it);

	cudaFree(d_t);cudaFree(d_V);cudaFree(d_dV2);cudaFree(d_Vnew);cudaFree(d_it);
	cudaFree(d_m);cudaFree(d_h);cudaFree(d_jj);cudaFree(d_d);
	cudaFree(d_f);cudaFree(d_X);cudaFree(d_cai);

	cudaFree(d_m0);cudaFree(d_h0);cudaFree(d_jj0);cudaFree(d_d0);
        cudaFree(d_f0);cudaFree(d_X0);cudaFree(d_dVdt);cudaFree(dcai);
}

void Send_to_Device(){
        cudaError_t Error;
        size_t size;
        size = nx*ny*sizeof(float);

	Error = cudaMemcpy(d_t, h_t, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_t->d_t) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(d_V, h_V, (nx+2)*(ny+2)*sizeof(float), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_V->d_V) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_m, h_m, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_m->d_m) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_h, h_h, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_h->d_h) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_jj, h_jj, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_jj->d_jj) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_d, h_d, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_d->d_d) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_f, h_f, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_f->d_f) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_X, h_X, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_X->d_X) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_cai, h_cai, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_cai->d_cai) = %s\n",cudaGetErrorString(Error));
}

void Send_to_Host(){
	cudaError_t Error;
        size_t size;
        size = nx*ny*sizeof(float);

        Error = cudaMemcpy(h_V, d_V, (nx+2)*(ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_V->h_V) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_m, d_m, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_m->h_m) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_h, d_h, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_h->h_h) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_jj, d_jj, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_jj->h_jj) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_d, d_d, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_d->h_d) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_f, d_f, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_f->h_f) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_X, d_X, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_X->h_X) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(h_cai, d_cai, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
	printf("CUDA error(copy d_cai->h_cai) = %s\n",cudaGetErrorString(Error));
}

void Send_V(){
        cudaError_t Error;
        size_t size;
        size = nx*ny*sizeof(float);

        Error = cudaMemcpy(h_Vnew, d_Vnew, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_Vnew->Vnew) = %s\n",cudaGetErrorString(Error));

}

void Save_Result(){

        FILE *pFile;
        int i,j;
        int index;
        //int n;
        //n = nx;
        pFile = fopen("V.txt","w+");
        // Save the matrix V
        for (i = 0; i < ny; i++) {
                for (j = 0; j < nx; j++) {
                        index = i*nx + j;
                        fprintf(pFile, "%g", h_Vnew[index]);
                        if (j == (nx-1)) {
                                fprintf(pFile, "\n");
                        }else{
                                fprintf(pFile, "\t");
                        }
                }
        }
        fclose(pFile);
/*
	pFile = fopen("Vnotnew.txt","w+");
        // Save the matrix V
        for (i = 1; i < ny+1; i++) {
                for (j = 1; j < nx+1; j++) {
                        index = i*(nx+2) + j;
                        fprintf(pFile, "%g", h_V[index]);
                        if (j == (nx)) {
                                fprintf(pFile, "\n");
                        }else{
                                fprintf(pFile, "\t");
                        }
                }
        }
        fclose(pFile);*/
}

