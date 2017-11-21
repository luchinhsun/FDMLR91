// BarModel1.cpp : 定義控制台應用程式的入口點。

#include <math.h>
//#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <sys/timeb.h>

//-*******FDM parameters for LR91 *******
//int const nx = 1000, ny = 1000;//grid numbers
//float dx = 0.015, dy = 0.015;//space step, 3cm*3cm
//float D = 0.001;//D: diffusion coefficient cm^2/ms
#define nx 200
#define ny 200
#define dx 0.015
#define dy 0.015
#define D 0.001

// Time Step
//float dt = 0.02; // Time step (ms)
//float t; // Time (ms)
//int steps; // Number of Steps
//int increment; // Loop Control Variable
//#define dt 0.02
#define dt_max 0.04
#define dt_min 0.001

// Voltage
//float V[nx + 2][nx + 2]; // Initial Voltage (mv)
//float dV2[nx + 2][nx + 2]; // second order derivatives of Voltage (mv)
//float Vnew[nx + 2][nx + 2];// New Voltage (mV)
//float dvdt; // Change in Voltage / Change in Time (mV/ms)
//float dvdtnew; // New dv/dt (mV/ms)

// Total Current and Stimulus
//float st; // Constant Stimulus (uA/cm^2)
//float tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
//float stimtime; //Time period during which stimulus is applied (ms)
//float it[nx + 1][nx + 1]; // Total current (uA/cm^2)
#define st (-80.0)
#define stimtime ((int)(0.6 / dt_max + 0.6))

// Terms for Solution of Conductance and Reversal Potential
//const float R = 8314; // Universal Gas Constant (J/kmol*K)
//const float frdy = 96485; // Faraday's Constant (C/mol)
//float temp = 310; // Temperature (K)
#define R 8314.0
#define frdy 96485.0
#define temp 310.0

// Ion Concentrations
//float nai; // Intracellular Na Concentration (mM)
//float nao; // Extracellular Na Concentration (mM)
//float cai[nx + 1][nx + 1]; // Intracellular Ca Concentration (mM)
//float cao; // Extracellular Ca Concentration (mM)
//float ki; // Intracellular K Concentration (mM)
//float ko; // Extracellular K Concentration (mM)
#define nai 18.0
#define nao 140.0
#define ki 145.0
#define ko 5.4
#define cao 1.8

/*
// Fast Sodium Current (time dependant)
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

// Current through L-type Ca Channel
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

// Time-dependent potassium current
float ik[nx + 1][nx + 1]; // Rapidly Activating K Current (uA/uF)
float gk; // Channel Conductance of Rapidly Activating K Current (mS/uF)
float ek; // Reversal Potential of Rapidly Activating K Current (mV)
float ax; // K alpha-x rate constant (ms^-1)
float bx; // K beta-x rate constant (ms^-1)
float X[nx + 1][nx + 1]; // Rapidly Activating K time-dependant activation  --gate X in LR91
float xss; // Steady-state value of inactivation gate xr  --gate X in LR91
float taux; // Time constant of gate xr (ms^-1) --gate X in LR91
float Xi; // K time-independent inactivation --gate Xi in LR91

// Potassium Current (time-independent)
float ik1[nx + 1][nx + 1]; // Time-independent K current (uA/uF)
float gk1; // Channel Conductance of Time Independant K Current (mS/uF)
float ek1; // Reversal Potential of Time Independant K Current (mV)
float ak1; // K alpha-ki rate constant (ms^-1)
float bk1; // K beta-ki rate constant (ms^-1)
float K1ss; // Steady-state value of K inactivation gate K1

// Plateau Potassium Current
float ikp[nx + 1][nx + 1]; // Plateau K current (uA/uF)
float gkp; // Channel Conductance of Plateau K Current (mS/uF)
float ekp; // Reversal Potential of Plateau K Current (mV)
float kp; // K plateau factor

// Background Current
float ib[nx + 1][nx + 1]; // Background current (uA/uF)

//performance compared
float Vmax, dvdt_max = 0, APD90, TNP, CPUtime;
float Vold, v_onset;
int flag = 0;
long f_count = 0;

// Ion Current Functions
void comp_ina(int i, int j); // Calculates Fast Na Current
void comp_ical(int i, int j); // Calculates Currents through L-Type Ca Channel
void comp_ik(int i, int j); // Calculates Time-dependent K Current
void comp_ik1(int i, int j); // Calculates Time-Independent K Current
void comp_ikp(int i, int j); // Calculates Plateau K Current
void comp_ib(int i, int j); // Calculates Background Current
void comp_it(int i, int j); // Calculates Total Current
*/

void Allocate();
void Send_to_Device();
void bc();
void dV2();
void dVdt();
void stimu();
void ODE_stim();
void ODE();
void Forward_Euler();
void Send_to_Host();
void Send_V();
void Save_Result();
void free();
