#pragma once
namespace myData
{
const double M_PI			= 3.141592653589793238462643;
//&&&&&&&&&&&&&&&&&

const double L				=	2*68.0;//Кинетичсекий момент маховика Hmc
const double wm_max			= 6000.0 * 2.0 * M_PI / 60.0;
const double M_max			= 170.0;
const double N_max			= 6000.0;

const double Jm             = L / wm_max; // момент инерции маховика kg/m^2

const double Jx				= 7986.343;//Проекции момента инерции KA
const double Jy				= 7987.281;//
const double Jz				= 3030.923;//

const double wx_min			= 0.001;//минимально допустимые угловые скорости
const double wy_min			= 0.001;//
const double wz_min			= 0.001;//

const double wx_start		= 0.35 * M_PI / 180;//Начальные угловые скорости КА
const double wy_start		= 0.7 * M_PI / 180;//
const double wz_start		= 0.7 * M_PI / 180;//

const double y_start		= 0.0 * M_PI / 180;//Начальное угловое положение KA
const double fi_start		= 0.0 * M_PI / 180;//
const double theta_start    = 0.0 * M_PI / 180;//

const double time_relax		=   100.0*60;//время успокоения, сek
//&&&&&&&&&&&&&&&&

const double w				= 2.0 * M_PI / time_relax;//собственная частота
const double e				=  sqrt(2) / 2.0;


const double J1				= Jx - Jm; //Момент инерции КА в отсутствие маховиков
const double J2				= Jy - Jm; //
const double J3				= Jz - Jm; //

const double Kpx			= J1 * w*w;
const double Kpy			= J2 * w*w;
const double Kpz			= J3 * w*w;

const double Krx			= 2 * J1*w*e;
const double Kry			= 2 * J2*w*e;
const double Krz			= 2 * J3*w*e;




const double Time_end = 14400.0;
}
#define ON     false
#define OFF    true
#define Radian OFF
#define First_Run 1