#pragma once

#define ON true
#define OFF false




#define light_Null 0
#define light_Midle 1
#define light_Full 2

#define MOON  ON
#define SUN   ON

#define First_run 1
#define Print ON


const bool print = ON;

const int Time_fly = (31 + 28 + 31) * 24; // в часах
const double step = 60.0;	// в секундах
const double eps = 100000.0;
const double R_GCO = 42164.0;
const double Time_start = 5.0 * 24.0 * 60.0 * 60.0;
const double RSun = 696000.0;//Sun radius
const double REarth = 6371.0;//Earth raus
const double Betta = 28.7*M_PI / 180;


const double Unit_R = R_GCO;//км
const double Unit_V = sqrt(nu_Earth / Unit_R);//км/с
const double Unit_T = Unit_R / Unit_V;//c
const double Unit_A = nu_Earth / (pow(Unit_R, 2));//км/с^2
const double nu_Earth1 = 1;

const double T_Earth = 86164.090530833;