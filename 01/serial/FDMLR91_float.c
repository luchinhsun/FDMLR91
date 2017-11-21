// BarModel1.cpp : 定義控制台應用程式的入口點。
//

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/timeb.h>

using std::cout;
using std::endl;

//-*******FDM parameters for LR91 *******
int const nx = 1000, ny = 1000;//grid numbers
float dx = 0.015, dy = 0.015;//space step, 3cm*3cm
float D = 0.001;//D: diffusion coefficient cm^2/ms

/* Time Step */
float dt = 0.02; // Time step (ms)
float t; // Time (ms)
int steps; // Number of Steps
int increment; // Loop Control Variable

/* Voltage */
float V[nx + 2][nx + 2]; // Initial Voltage (mv)
float dV2[nx + 2][nx + 2]; // second order derivatives of Voltage (mv)
float Vnew[nx + 2][nx + 2];// New Voltage (mV)
float dvdt; // Change in Voltage / Change in Time (mV/ms)
float dvdtnew; // New dv/dt (mV/ms)

/* Total Current and Stimulus */
float st; // Constant Stimulus (uA/cm^2)
float tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
float stimtime; //Time period during which stimulus is applied (ms)
float it[nx + 1][nx + 1]; // Total current (uA/cm^2)

/* Terms for Solution of Conductance and Reversal Potential */
const float R = 8314; // Universal Gas Constant (J/kmol*K)
const float frdy = 96485; // Faraday's Constant (C/mol)
float temp = 310; // Temperature (K)

/* Ion Concentrations */
float nai; // Intracellular Na Concentration (mM)
float nao; // Extracellular Na Concentration (mM)
float cai[nx + 1][nx + 1]; // Intracellular Ca Concentration (mM)
float cao; // Extracellular Ca Concentration (mM)
float ki; // Intracellular K Concentration (mM)
float ko; // Extracellular K Concentration (mM)

/* Fast Sodium Current (time dependant) */
float ina[nx + 1][nx + 1]; // Fast Na Current (uA/uF)
float gna; // Max. Conductance of the Na Channel (mS/uF)
float ena; // Reversal Potential of Na (mV)
float am; // Na alpha-m rate constant (ms^-1)
float bm; // Na beta-m rate constant (ms^-1)
float ah; // Na alpha-h rate constant (ms^-1)
float bh; // Na beta-h rate constant (ms^-1)
float aj; // Na alpha-j rate constant (ms^-1)
float bj; // Na beta-j rate constant (ms^-1)
float mtau; // Na activation
float htau; // Na inactivation
float jtau; // Na inactivation
float mss; // Na activation
float hss; // Na inactivation
float jss; // Na slow inactivation
float m[nx + 1][nx + 1]; // Na activation
float h[nx + 1][nx + 1]; // Na inactivation
float jj[nx + 1][nx + 1]; // Na slow inactivation

/* Current through L-type Ca Channel */
float dcai; // Change in myoplasmic Ca concentration (mM)
float isi[nx + 1][nx + 1]; // Slow inward current (uA/uF)
float esi[nx + 1][nx + 1]; // Reversal Potential of si (mV)
float ad; // Ca alpha-d rate constant (ms^-1)
float bd; // Ca beta-d rate constant (ms^-1)
float af; // Ca alpha-f rate constant (ms^-1)
float bf; // Ca beta-f rate constant (ms^-1)

float d[nx + 1][nx + 1]; // Voltage dependant activation gate
float dss; // Steady-state value of activation gate d
float taud; // Time constant of gate d (ms^-1)----mistake ???？ms？
float f[nx + 1][nx + 1]; // Voltage dependant inactivation gate
float fss; // Steady-state value of inactivation gate f
float tauf; // Time constant of gate f (ms^-1)
float fca[nx + 1][nx + 1]; // Ca dependant inactivation gate -from LR94

/* Time-dependent potassium current*/
float ik[nx + 1][nx + 1]; // Rapidly Activating K Current (uA/uF)
float gk; // Channel Conductance of Rapidly Activating K Current (mS/uF)
float ek; // Reversal Potential of Rapidly Activating K Current (mV)
float ax; // K alpha-x rate constant (ms^-1)
float bx; // K beta-x rate constant (ms^-1)
float X[nx + 1][nx + 1]; // Rapidly Activating K time-dependant activation  --gate X in LR91
float xss; // Steady-state value of inactivation gate xr  --gate X in LR91
float taux; // Time constant of gate xr (ms^-1) --gate X in LR91
float Xi; // K time-independent inactivation --gate Xi in LR91

/* Potassium Current (time-independent) */
float ik1[nx + 1][nx + 1]; // Time-independent K current (uA/uF)
float gk1; // Channel Conductance of Time Independant K Current (mS/uF)
float ek1; // Reversal Potential of Time Independant K Current (mV)
float ak1; // K alpha-ki rate constant (ms^-1)
float bk1; // K beta-ki rate constant (ms^-1)
float K1ss; // Steady-state value of K inactivation gate K1

/* Plateau Potassium Current */
float ikp[nx + 1][nx + 1]; // Plateau K current (uA/uF)
float gkp; // Channel Conductance of Plateau K Current (mS/uF)
float ekp; // Reversal Potential of Plateau K Current (mV)
float kp; // K plateau factor

/* Background Current */
float ib[nx + 1][nx + 1]; // Background current (uA/uF)

//performance compared
float Vmax, dvdt_max = 0, APD90, TNP, CPUtime;
float Vold, v_onset;
int flag = 0;
long f_count = 0;

/* Ion Current Functions */
void comp_ina(int i, int j); // Calculates Fast Na Current
void comp_ical(int i, int j); // Calculates Currents through L-Type Ca Channel
void comp_ik(int i, int j); // Calculates Time-dependent K Current
void comp_ik1(int i, int j); // Calculates Time-Independent K Current
void comp_ikp(int i, int j); // Calculates Plateau K Current
void comp_ib(int i, int j); // Calculates Background Current
void comp_it(int i, int j); // Calculates Total Current

int main(int argc, char* argv[])
{
	/* Data File */
	FILE *ap;
	FILE *fevaluation;
	fevaluation = fopen("fevaluation", "w");

	/* Time Loop Conditions */
	t = 0.0; // Time (ms)
	//	steps = (bcl*beats)/udt; // Number of ms
	st = -80.0; // Stimulus (mA)

	/* Beginning Ion Concentrations */
	nai = 18; // Initial Intracellular Na (mM)
	nao = 140; // Initial Extracellular Na (mM)
	ki = 145; // Initial Intracellular K (mM)
	ko = 5.4; // Initial Extracellular K (mM)
	//cai = 0.0002; // Initial Intracellular Ca (mM)
	cao = 1.8; // Initial Extracellular Ca (mM)

	int ncount, i, j;

	for (i = 1; i < nx + 1; i++){
		for (j = 1; j < ny + 1; j++){
			V[i][j] = -88.654973; // Initial Voltage (mv)
			m[i][j] = 0.000838;
			h[i][j] = 0.993336;
			jj[i][j] = 0.995484;
			d[i][j] = 0.000003;
			f[i][j] = 0.999745;
			X[i][j] = 0.000129;
			cai[i][j] = 0.0002; // Initial Intracellular Ca (mM)
		}
	}

	int nstep = 500; // snapshot interval to save data files 500*0.02=10 ms
	int index = 0;// filename index from 1-5
	char filename[100];

	struct timeb start, end;
        int diff;
        ftime(&start);
	for (ncount = 0; ncount <= 50000; ncount++){//30000 steps, 600ms
		for (i = 1; i < nx + 1; i++){
			//-****no flux boundary conditions*****
			V[i][0] = V[i][1];
			V[i][ny + 1] = V[i][ny];
			for (j = 1; j < ny + 1; j++){
				V[0][j] = V[1][j];
				V[nx + 1][j] = V[nx][j];
			}
		}
		//-*********** Center Differnce for Space *******
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				comp_ina(i, j);
				comp_ical(i, j);
				comp_ik(i, j);
				comp_ik1(i, j);
				comp_ikp(i, j);
				comp_ib(i, j);
				comp_it(i, j);

				dV2[i][j] = -it[i][j] + D*((V[i + 1][j] + V[i - 1][j] - 2 * V[i][j]) / (dx*dx) +
						(V[i][j + 1] + V[i][j - 1] - 2 * V[i][j]) / (dy*dy));
			}
		}

		//-*****stimulation with a plane waves****
		if (ncount >= 0 && ncount <= 100) { //stimulus is hold with 0.5 ms, 0.02*25 = 0.5ms
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j <= 5; j++){//最少3列細胞才能激發平面波
					dV2[i][j] = dV2[i][j] + (-st);
				}
			}
		}

		int fileflag = 0;

		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				//Forward Euler
				Vnew[i][j] = V[i][j] + dt*dV2[i][j];
				V[i][j] = Vnew[i][j];
				if (ncount%nstep == 0){//get data at the 9000th step
					if (fileflag == 0){
						sprintf(filename, "ap%d", index);
						ap = fopen(filename, "w");
						fileflag = 1;
						index++;
					}
					fprintf(ap, "%g\t", V[i][j]);
					if (j == ny){
						fprintf(ap, "\n");
					}
				}
			}
		}
		if (fileflag == 1){
			fclose(ap);
		}

		t = t + dt;
		//-***********trancation 1/2 of the plane wave to generate a spiral wave******
		if (ncount == 20000){
			for (i = 1; i < nx / 2; i++){
				for (j = 1; j < ny; j++){
					V[i][j] = -88.654973; // Initial Voltage (mv)
					m[i][j] = 0.000838;
					h[i][j] = 0.993336;
					jj[i][j] = 0.995484;
					d[i][j] = 0.000003;
					f[i][j] = 0.999745;
					X[i][j] = 0.000129;
					cai[i][j] = 0.0002; // Initial Intracellular Ca (mM)
				}
			}
		}
	}
	ftime(&end);
        diff = (int)(1000.0*(end.time-start.time)+(end.millitm-start.millitm));
        printf("\nTime = %d ms\n", diff);
	//float time_used = (float)(end - start) / CLK_TCK;
	//fprintf(fevaluation, "%g", time_used);
}
/********************************************************/
/* Functions that describe the currents begin here */

//Fast sodium current
void comp_ina(int i, int j) {
	gna = 23;
	ena = ((R*temp) / frdy)*log(nao / nai);

	am = 0.32*(V[i][j] + 47.13) / (1 - exp(-0.1*(V[i][j] + 47.13)));
	bm = 0.08*exp(-V[i][j] / 11);
	if (V[i][j] < -40) {
		ah = 0.135*exp((80 + V[i][j]) / -6.8);
		bh = 3.56*exp(0.079*V[i][j]) + 310000 * exp(0.35*V[i][j]);
		aj = (-127140 * exp(0.2444*V[i][j]) - 0.00003474*exp(-0.04391*V[i][j]))*((V[i][j] + 37.78) / (1 + exp(0.311*(V[i][j] + 79.23))));
		bj = (0.1212*exp(-0.01052*V[i][j])) / (1 + exp(-0.1378*(V[i][j] + 40.14)));
	}
	else {
		ah = 0;
		bh = 1 / (0.13*(1 + exp((V[i][j] + 10.66) / -11.1)));
		aj = 0;
		bj = (0.3*exp(-0.0000002535*V[i][j])) / (1 + exp(-0.1*(V[i][j] + 32)));
	}
	mtau = 1 / (am + bm);
	htau = 1 / (ah + bh);
	jtau = 1 / (aj + bj);

	mss = am*mtau;
	hss = ah*htau;
	jss = aj*jtau;

	m[i][j] = mss - (mss - m[i][j])*exp(-dt / mtau);
	h[i][j] = hss - (hss - h[i][j])*exp(-dt / htau);
	jj[i][j] = jss - (jss - jj[i][j])*exp(-dt / jtau);

	ina[i][j] = gna*m[i][j] * m[i][j] * m[i][j] * h[i][j] * jj[i][j] * (V[i][j] - ena);
}

//Slow inward current
void comp_ical(int i, int j) {
	esi[i][j] = 7.7 - 13.0287*log(cai[i][j]);

	ad = 0.095*exp(-0.01*(V[i][j] - 5)) / (1 + exp(-0.072*(V[i][j] - 5)));
	bd = 0.07*exp(-0.017*(V[i][j] + 44)) / (1 + exp(0.05*(V[i][j] + 44)));
	af = 0.012*exp(-0.008*(V[i][j] + 28)) / (1 + exp(0.15*(V[i][j] + 28)));
	bf = 0.0065*exp(-0.02*(V[i][j] + 30)) / (1 + exp(-0.2*(V[i][j] + 30)));

	taud = 1 / (ad + bd);
	tauf = 1 / (af + bf);

	dss = ad*taud;
	fss = af*tauf;

	d[i][j] = dss - (dss - d[i][j])*exp(-dt / taud);
	f[i][j] = fss - (fss - f[i][j])*exp(-dt / tauf);

	isi[i][j] = 0.09*d[i][j] * f[i][j] * (V[i][j] - esi[i][j]);

	dcai = -0.0001*isi[i][j] + 0.07*(0.0001 - cai[i][j]);

	cai[i][j] = cai[i][j] + dcai*dt;//Ca的變化量
}

//Time-dependent potassium current
void comp_ik(int i, int j) {
	gk = 0.282*sqrt(ko / 5.4);
	ek = ((R*temp) / frdy)*log(ko / ki);
	//float prnak = 0.01833;
	//ek = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));

	ax = 0.0005*exp(0.083*(V[i][j] + 50)) / (1 + exp(0.057*(V[i][j] + 50)));
	bx = 0.0013*exp(-0.06*(V[i][j] + 20)) / (1 + exp(-0.04*(V[i][j] + 20)));

	taux = 1 / (ax + bx);
	xss = ax*taux;
	X[i][j] = xss - (xss - X[i][j])*exp(-dt / taux);

	if (V[i][j] > -100) {
		Xi = 2.837*(exp(0.04*(V[i][j] + 77)) - 1) / ((V[i][j] + 77 + 1e-15)*exp(0.04*(V[i][j] + 35)));
	}
	else {
		Xi = 1;
	}

	ik[i][j] = gk*X[i][j] * Xi*(V[i][j] - ek);
}


//Time-independent potassium current
void comp_ik1(int i, int j) {
	gk1 = 0.6047*(sqrt(ko / 5.4));
	ek1 = ((R*temp) / frdy)*log(ko / ki);

	ak1 = 1.02 / (1 + exp(0.2385*(V[i][j] - ek1 - 59.215)));
	bk1 = (0.49124*exp(0.08032*(V[i][j] - ek1 + 5.476)) + exp(0.06175*(V[i][j] - ek1 - 594.31)))
			/(1 + exp(-0.5143*(V[i][j] - ek1 + 4.753)));
	K1ss = ak1 / (ak1 + bk1);

	ik1[i][j] = gk1*K1ss*(V[i][j] - ek1);
}

//Plateau potassium current
void comp_ikp(int i, int j) {
	gkp = 0.0183;
	ekp = ek1;

	kp = 1 / (1 + exp((7.488 - V[i][j]) / 5.98));

	ikp[i][j] = gkp*kp*(V[i][j] - ekp);
}

//Background current
void comp_ib(int i, int j) {
	ib[i][j] = 0.03921*(V[i][j] + 59.87);
}

/* Total sum of currents is calculated here, if the time is between
stimtime = 0 and stimtime = 0.5 (ms), a stimulus is applied */
//%刺激電流的持續時間限制在0-0.5之間，超過刺激電流就置零
void comp_it(int i, int j) {
	//當時間t到達10.01ms後，刺激電流才引入
	//
	//	if (t >= 5 && t<(5 + 0.5)) {
	//		it[i][j] = st + ina[i][j] + isi[i][j] + ik[i][j] + ik1[i][j] + ikp[i][j] + ib[i][j];
	//	}else {
	it[i][j] = ina[i][j] + isi[i][j] + ik[i][j] + ik1[i][j] + ikp[i][j] + ib[i][j];
	//	}
}


/* Values are printed to a file called ap. The voltage and
currents can be plotted versus time using graphing software. */
//void prttofile() {
//	if (t>(0) && t<(bcl*beats))
//	{
//		fprintf(ap, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
//			t, v, nai, ki, cai, ina, isi, ikr, iki, ikp, ib);
//		//printf("%.5f\t%g\n", t, v);
//		//printf("%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
//		//	t, v, nai, ki, cai, ina, isi, ikr, iki, ikp, ib);
//	}
//	//nai, ki, cai are the Intracellular Concentration of nai, ki, cai
//}

