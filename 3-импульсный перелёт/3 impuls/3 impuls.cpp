
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <vector>
#include <thread>

#include "Planet_position.h"
#include "ConstData.h"


class Timer
{
public:
	Timer()
	{
		start = std::chrono::high_resolution_clock::now();
	}
	~Timer()
	{
		end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> difer = end - start;
		std::cout << " Duratioun - " << difer.count() << std::endl;
	}
private:
	std::chrono::time_point<std::chrono::steady_clock> start, end;
};




class KA
{
public:
	static int count;

public:
	double  x, y, z;
	double Vx, Vy, Vz, VLx4, VLy4, VLz4;
	double x2, y2, z2, Vx2, Vy2, Vz2, dVx2, dVy2, dVz2;
	double xL, yL, zL;
	double delV, delV3, delV2;
	double Time, d, Flight_Time1, Time_start, Time_end, Ful_time, T2, OptimalT;
	double V;
	int metka, point;
	double i, w, p, omeg, e, u;
	double Anom;

	KA difer(double Time_end1);
	KA Back_difer();
	KA Start_Position();
	int MMM;
	double Betaconst;	double lamdaconst;	 double tconst;

	double  MainOptimization(double dV3);
	void	OptimizationSpeed();
	double Radius();
	double  FunFlight(double Fly_t);

	double   Rynge_Kyt();
	double   Full_Rynge_Kyt();
	void	 Back_Rynge_Kyt();
	void     Rynge_Kytta_L4();
	void	 Degree(std::ofstream *fx, std::ofstream *fy, std::ofstream *fz);
	void	 print(std::ofstream *f, double Time_start);
	void	 Grafic_X(std::ofstream *f);
	void	 Grafic_Y(std::ofstream *f);
	void	 Grafic_Z(std::ofstream *f);
	void	 Grafic_Orbit_Parametr(std::ofstream *e, std::ofstream *i, std::ofstream *omeg, std::ofstream *w, std::ofstream *a, std::ofstream *h);
	void	 Track(std::ofstream *lam, std::ofstream *bet);
	void	 ShowParametr();
	double   E();


	double Y1, W2, m_adap, m_rb_D, m_ka_after_first_ipmusl, m_ka_after_second_impuls, m_fuil_for_first_impuls, m_fuil_for_second_ipmuls, m_NOO, m_1, m_2, m_3;

};

int KA::count = 0;
#include "Input.h"
double KA::Radius()
{
	return sqrt(x * x + y * y + z * z);
}






void Mult_Matrix(double **a, double *b, int size)
{
	double *x = new double[size];
	double sum = 0;
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
			sum += a[j][i] * b[i];
		x[j] = sum;
		sum = 0;
	}

	for (int i = 0; i < size; i++)
		b[i] = x[i];
	delete[] x;
}

void KA::Degree(std::ofstream *fx, std::ofstream *fy, std::ofstream *fz)
{

	int j;
	double hi = Time * Unit_T;

	for (size_t i = 0; i < Month; i++)
	{
		if (hi / (60) < H[i + 1] && hi / (60.0) > H[i])
		{
			j = i;
			break;
		}

		else if (hi / (60) == H[i])
		{
			j = i;
			break;
		}

	}



	double Xl	= ((xL4[j] + (hi / (60) - H[j])*(xL4[j + 1] - xL4[j]) / (H[j + 1] - H[j])));
	double Yl	= ((yL4[j] + (hi / (60) - H[j])*(yL4[j + 1] - yL4[j]) / (H[j + 1] - H[j])));
	double Zl	= ((zL4[j] + (hi / (60) - H[j])*(zL4[j + 1] - zL4[j]) / (H[j + 1] - H[j])));
	double mod  = sqrt(pow(Xl, 2) + pow(Yl, 2) + pow(Zl, 2));
	double deg;
	double l;

	l = asin((Zl) / (mod));


	//std::cout <<" Time = "<<Time* Unit_T<<" X = "<<Xl<< " Y = "<<Yl<<" Z = "<<Zl<<" Betta = " << asin(Zl / mod)*180/M_PI << std::endl;
	if (Yl >= 0)
		deg = (acos((Xl) / (mod*cos(l))));// *cos(Betta))))*180.0 / M_PI;
	else
		deg = 2.0 *M_PI - (acos((Xl) / (mod*cos(l))));// *cos(Betta))))*180.0 / M_PI;

	//std::cout << " deg = " << deg*180/M_PI << " Betta = " << asin(Zl / mod) * 180 / M_PI << std::endl;
	double **A = new double *[3];
	for (int i = 0; i < 3; i++)
		A[i] = new double[3];

	double *Vec = new double[3];

	A[0][0] = cos(l)*cos(deg);		A[0][1] = cos(l)*sin(deg);		A[0][2] = sin(l);
	A[1][0] = -sin(deg);			A[1][1] = cos(deg);				A[1][2] = 0;
	A[2][0] = -cos(deg)*sin(l);		A[2][1] = -sin(l)*sin(deg);		A[2][2] = cos(l);

	Vec[0] = x;	Vec[1] = y;	Vec[2] = z;

	Mult_Matrix(A, Vec, 3);


	*fx << Vec[0] * Unit_R << ",";	*fy << Vec[1] * Unit_R << ",";	*fz << Vec[2] * Unit_R << ",";


	for (int i = 0; i < 3; i++)
		delete A[i];
	delete[] A;

	delete[] Vec;


}

void KA::Track(std::ofstream *lam, std::ofstream *bet)
{
	double Gx	= y * Unit_R *Vz*Unit_V - z * Unit_R*Vy*Unit_V;
	double Gy	= z * Unit_R*Vx*Unit_V - x * Unit_R*Vz*Unit_V;
	double Gz	= x * Unit_R*Vy*Unit_V - y * Unit_R*Vx*Unit_V;

	double G	= sqrt(Gx*Gx + Gy * Gy + Gz * Gz);
	double v	= sqrt(pow(Vx*Unit_V, 2) + pow(Vy*Unit_V, 2) + pow(Vz*Unit_V, 2));
	double r	= sqrt(pow(x*Unit_R, 2) + pow(y*Unit_R, 2) + pow(z*Unit_R, 2));

	double h	= v * v - (2 * nu_Earth / r);
	double e	= sqrt(1 + ((G*G)*h / (nu_Earth*nu_Earth)));

	double Gz0	= (Gz / G);
	double Gy0	= (Gy / G);
	double Gx0	= (Gx / G);

	double i = acos(Gz0);

	double omeg;
	if (Gx0 >= 0)
		omeg = acos((-Gy0) / sin(i));
	else
		omeg = 2 * M_PI - acos((-Gy0) / sin(i));

	double lx = (-(Gy*Vz*Unit_V - Gz * Vy*Unit_V)) - (nu_Earth*x*Unit_R / r);
	double ly = (-(Gz*Vx*Unit_V - Gx * Vz*Unit_V)) - (nu_Earth*y*Unit_R / r);
	double lz = (-(Gx*Vy*Unit_V - Gy * Vx*Unit_V)) - (nu_Earth*z*Unit_R / r);

	double l = sqrt(lx*lx + ly * ly + lz * lz);

	double w;

	if ((lz / l) > 0)
	{
		w = acos(((lx / l)*cos(omeg)) + ((ly / l)*sin(omeg)));
	}
	else
	{
		w = (2 * M_PI) - acos(((lx / l)*cos(omeg)) + ((ly / l)*sin(omeg)));
	}

	double p	= (G*G) / nu_Earth;
	double a	= p / (1 - pow(e, 2));

	double ckolar;
	double anom, x1;
	ckolar = x * Vx + y * Vy + z * Vz;

	if (ckolar >= 0)
	{

		x1 = ((p / r) - 1) / e;

		anom = acos(x1);
	}
	else
	{
		anom = 2 * M_PI - acos(((p / r) - 1) / e);
	}


	double Betta;	double lamda;	double fi;
	double wr = 2 * M_PI / T_Earth;

	fi = asin(z / Radius());
	if (x > 0)
		Betta = atan(y / x);
	else if (x < 0 && y >= 0)
		Betta = atan(y / x) + M_PI;
	else if (x < 0 && y < 0)
		Betta = atan(y / x) - M_PI;
	else if (x = 0 && y > 0)
		Betta = M_PI / 2;
	else if (x = 0 && y < 0)
		Betta = -(M_PI / 2);
	else
		std::cout << "мы над полюсом" << std::endl;

	//lamda = - wr*
	//*lam<<
	//lamda = Betta-wr*
}

double KA::E()
{

	double *min = new double[6];

	min[0] = x;		min[1] = Vx;
	min[2] = y;		min[3] = Vy;
	min[4] = z;		min[5] = Vz;
	/*
	if (x > y && x < z)
		return x;
	else if (y > x && y > z)
		return y;
	else return z;*/
	double max = min[0];
	for (int i = 1; i < 6; i++)
	{
		if (min[i] > max)
			max = min[i];
	}


	delete[] min;
	return max;
}

KA KA::difer(double Time_end1)
{
	count++;









	double Mx; double My; double Mz;
	double Sx, Sy, Sz;


	int j;
	double hi = Time * Unit_T;
	int End = Time_end1/(60*60.0) + 20;
	for (size_t i = 0; i < End; i++)
	{
		if (hi / (60) < H[i + 1] && hi / (60) > H[i])
		{
			j = i;
			break;
		}

		else if (hi / (60) == H[i])
		{
			j = i;
			break;

		}

	}


	//std::cout << " j = " << j << " H = " << H[j] << std::endl;
	Mx = (xM[j] + (hi / (60) - H[j])*(xM[j + 1] - xM[j]) / (H[j + 1] - H[j])) / R_GCO;
	My = (yM[j] + (hi / (60) - H[j])*(yM[j + 1] - yM[j]) / (H[j + 1] - H[j])) / R_GCO;
	Mz = (zM[j] + (hi / (60) - H[j])*(zM[j + 1] - zM[j]) / (H[j + 1] - H[j])) / R_GCO;

	Sx = (xS[j] + (hi / (60) - H[j])*(xS[j + 1] - xS[j]) / (H[j + 1] - H[j])) / R_GCO;
	Sy = (yS[j] + (hi / (60) - H[j])*(yS[j + 1] - yS[j]) / (H[j + 1] - H[j])) / R_GCO;
	Sz = (zS[j] + (hi / (60) - H[j])*(zS[j + 1] - zS[j]) / (H[j + 1] - H[j])) / R_GCO;



	double t = Time;


	double R_KA = (sqrt(x*x + y * y + z * z));
	double R_KA_M = ((sqrt((x - Mx)*(x - Mx) + (y - My)*(y - My) + (z - Mz)*(z - Mz))));
	double R_M = (sqrt(Mx*Mx + My * My + Mz * Mz));

	double R_KA_S = ((sqrt((x - Sx)*(x - Sx) + (y - Sy)*(y - Sy) + (z - Sz)*(z - Sz))));
	double R_S = (sqrt(Sx*Sx + Sy * Sy + Sz * Sz));

	double nu_Moon1 = nu_Moon / nu_Earth;
	double nu_Sun1 = nu_Sun / nu_Earth;


	//double Unit_a = nu_Earth / 42164.0;

	KA d;

	d.x = Vx;
	d.y = Vy;
	d.z = Vz;




	if (SUN && MOON)
	{
		d.Vx = (-nu_Earth1 * x / (pow(R_KA, 3))) - (nu_Moon1 * (x - Mx) / (pow(R_KA_M, 3))) - (nu_Moon1 * Mx / (pow(R_M, 3))) - (nu_Sun1 * (x - Sx) / (pow(R_KA_S, 3))) - (nu_Sun1 * Sx / (pow(R_S, 3)));// / Unit_a;
		d.Vy = (-nu_Earth1 * y / (pow(R_KA, 3)) - nu_Moon1 * (y - My) / (pow(R_KA_M, 3)) - nu_Moon1 * My / (pow(R_M, 3)) - nu_Sun1 * (y - Sy) / (pow(R_KA_S, 3)) - nu_Sun1 * Sy / (pow(R_S, 3)));
		d.Vz = -nu_Earth1 * z / (pow(R_KA, 3)) - nu_Moon1 * (z - Mz) / (pow(R_KA_M, 3)) - nu_Moon1 * Mz / (pow(R_M, 3)) - nu_Sun1 * (z - Sz) / (pow(R_KA_S, 3)) - nu_Sun1 * Sz / (pow(R_S, 3));
	}
	else if (MOON)
	{
		d.Vx = -nu_Earth * x / (pow(R_KA, 3)) - nu_Moon * (x - Mx) / (pow(R_KA_M, 3)) - nu_Moon * Mx / (pow(R_M, 3));
		d.Vy = -nu_Earth * y / (pow(R_KA, 3)) - nu_Moon * (y - My) / (pow(R_KA_M, 3)) - nu_Moon * My / (pow(R_M, 3));
		d.Vz = -nu_Earth * z / (pow(R_KA, 3)) - nu_Moon * (z - Mz) / (pow(R_KA_M, 3)) - nu_Moon * Mz / (pow(R_M, 3));
	}
	else if (MOON == OFF && SUN == OFF)
	{
		d.Vx = -nu_Earth1 * x / (pow(R_KA, 3));
		d.Vy = -nu_Earth1 * y / (pow(R_KA, 3));
		d.Vz = -nu_Earth1 * z / (pow(R_KA, 3));
	}



	d.Time = 1;

	if (Time*Unit_T == Time_fly * 60 * 60)
	{
		count = 0;
		j = 0;

	}

	return d;
}

KA KA::Back_difer()
{





	double Mx; double My; double Mz;
	double Sx, Sy, Sz;


	int j;
	double hi = Time * Unit_T;


	for (size_t i = 0; i < Month; i++)
	{
		if (hi / (60) < H[i + 1] && hi / (60) > H[i])
			j = i;
		else if (hi / (60) == H[i])
			j = i;
	}


	//std::cout << " j = " << j << " H = " << H[j] << std::endl;
	Mx = (xM[j] + (hi / (60) - H[j])*(xM[j + 1] - xM[j]) / (H[j + 1] - H[j])) / R_GCO;
	My = (yM[j] + (hi / (60) - H[j])*(yM[j + 1] - yM[j]) / (H[j + 1] - H[j])) / R_GCO;
	Mz = (zM[j] + (hi / (60) - H[j])*(zM[j + 1] - zM[j]) / (H[j + 1] - H[j])) / R_GCO;

	Sx = (xS[j] + (hi / (60) - H[j])*(xS[j + 1] - xS[j]) / (H[j + 1] - H[j])) / R_GCO;
	Sy = (yS[j] + (hi / (60) - H[j])*(yS[j + 1] - yS[j]) / (H[j + 1] - H[j])) / R_GCO;
	Sz = (zS[j] + (hi / (60) - H[j])*(zS[j + 1] - zS[j]) / (H[j + 1] - H[j])) / R_GCO;



	double t = Time;


	double R_KA = (sqrt(x*x + y * y + z * z));
	double R_KA_M = ((sqrt((x - Mx)*(x - Mx) + (y - My)*(y - My) + (z - Mz)*(z - Mz))));
	double R_M = (sqrt(Mx*Mx + My * My + Mz * Mz));

	double R_KA_S = ((sqrt((x - Sx)*(x - Sx) + (y - Sy)*(y - Sy) + (z - Sz)*(z - Sz))));
	double R_S = (sqrt(Sx*Sx + Sy * Sy + Sz * Sz));

	double nu_Moon1 = nu_Moon / nu_Earth;
	double nu_Sun1 = nu_Sun / nu_Earth;


	//double Unit_a = nu_Earth / 42164.0;

	KA d;

	d.x = -Vx;
	d.y = -Vy;
	d.z = -Vz;



	if (SUN && MOON)
	{
		d.Vx = -((-nu_Earth1 * x / (pow(R_KA, 3))) - (nu_Moon1 * (x - Mx) / (pow(R_KA_M, 3))) - (nu_Moon1 * Mx / (pow(R_M, 3))) - (nu_Sun1 * (x - Sx) / (pow(R_KA_S, 3))) - (nu_Sun1 * Sx / (pow(R_S, 3))));// / Unit_a;
		d.Vy = -((-nu_Earth1 * y / (pow(R_KA, 3)) - nu_Moon1 * (y - My) / (pow(R_KA_M, 3)) - nu_Moon1 * My / (pow(R_M, 3)) - nu_Sun1 * (y - Sy) / (pow(R_KA_S, 3)) - nu_Sun1 * Sy / (pow(R_S, 3))));
		d.Vz = -(-nu_Earth1 * z / (pow(R_KA, 3)) - nu_Moon1 * (z - Mz) / (pow(R_KA_M, 3)) - nu_Moon1 * Mz / (pow(R_M, 3)) - nu_Sun1 * (z - Sz) / (pow(R_KA_S, 3)) - nu_Sun1 * Sz / (pow(R_S, 3)));
	}
	else if (MOON)
	{
		d.Vx = (-nu_Earth * x / (pow(R_KA, 3)) - nu_Moon * (x - Mx) / (pow(R_KA_M, 3)) - nu_Moon * Mx / (pow(R_M, 3)));
		d.Vy = (-nu_Earth * y / (pow(R_KA, 3)) - nu_Moon * (y - My) / (pow(R_KA_M, 3)) - nu_Moon * My / (pow(R_M, 3)));
		d.Vz = -(-nu_Earth * z / (pow(R_KA, 3)) - nu_Moon * (z - Mz) / (pow(R_KA_M, 3)) - nu_Moon * Mz / (pow(R_M, 3)));
	}
	else if (MOON == OFF && SUN == OFF)
	{
		d.Vx = -(-nu_Earth1 * x / (pow(R_KA, 3)));
		d.Vy = -(-nu_Earth1 * y / (pow(R_KA, 3)));
		d.Vz = -(-nu_Earth1 * z / (pow(R_KA, 3)));
	}



	d.Time = -1;

	if (Time*Unit_T == Time_fly * 60 * 60)
	{
		count = 0;
		j = 0;

	}

	return d;
}

KA KA::Start_Position()
{
	double Vr; double Vn;	//радиальная и нормальная составляющая скорости
	double R;

	KA a1;

	Anom = u - w;
	std::cout << " R = " << p / (1 + e * cos(Anom));

	Vr = sqrt(nu_Earth / p)*e*sin(Anom);
	Vn = sqrt(nu_Earth / p)*(1 + e * cos(Anom));
	std::cout << " Vr= " << Vr << " Vтрансвер= " << Vn << std::endl;
	R = p / (1 + e * cos(Anom));
	a1.x = R * (cos(omeg)*cos(u) - sin(omeg)*sin(u)*cos(i));
	a1.y = R * (sin(omeg)*cos(u) + cos(omeg)*sin(u)*cos(i));
	a1.z = R * sin(i)*sin(u);

	a1.Vx = Vr * (cos(u)*cos(omeg) - sin(u)*cos(i)*sin(omeg)) + Vn * (-sin(u)*cos(omeg) - cos(i)*sin(omeg)*cos(u));
	a1.Vy = Vr * (sin(omeg) * cos(u) + cos(omeg) * sin(u) * cos(i)) + Vn * (-sin(u)*sin(omeg) + cos(u)*cos(i)*cos(omeg));
	a1.Vz = Vr * sin(i)*sin(u) + Vn * cos(u)*sin(i);

	return a1;
}

double Kvadrat(double x, double y, double z)
{
	return sqrt(x*x + y * y + z * z);
}


double Lighting(double x, double y, double z, double XSun, double YSun, double ZSun)
{
	double R_ka, rst, rst2, alfa, fi, beta, eta, rs, penumbra, umbra, delta, gamma;

	double Rnt	= REarth;
	double Rsun = RSun;

	rs		= Kvadrat(XSun, YSun, ZSun);
	R_ka	= Kvadrat(x, y, z);

	alfa		= acos((x*XSun + y * YSun + z * ZSun) / (R_ka*rs));
	eta			= acos((Rnt + Rsun) / rs);
	delta		= alfa - eta;
	umbra		= 0;
	penumbra	= 0;

	if (delta == M_PI / 2.0)
	{
		penumbra = 0.0;
		//std::cout << " sclmdc\n";

	}
	else
	{
		if (delta < M_PI / 2.0 && delta >= 0)
		{
			rst2 = Rnt / cos(delta);
			penumbra = (R_ka - rst2)*cos(delta);
		}
		else
		{
			rst2 = Rnt / cos(delta);
			penumbra = fabs((R_ka - rst2)*cos(delta));
		}

		if (delta > M_PI / 2.0 && delta <= M_PI)
		{
			rst2		= Rnt / cos(delta);
			penumbra	= -fabs((R_ka - rst2)*cos(delta));
		}
	}




	gamma		= acos((Rsun - Rnt) / rs);
	fi			= M_PI - gamma;
	beta		= alfa - fi;

	if (beta >= 0){
		rst		= Rnt / acos(beta);
		umbra	= (R_ka - rst)*cos(beta);
	}
	else{
		rst		= Rnt / acos(beta);
		umbra   = fabs((R_ka - rst)*cos(beta));
	}
	//std::cout << " penubmra = " << penumbra << " umbra = " << umbra << std::endl;

	if (umbra < 0)
		return 0;
	else if (penumbra < 0)
		return 1;
	else
		return 2;
}







double KA::Rynge_Kyt() //интегрирование до момента приложения второго импульса 
{

	double Vr; double Vn; double R;	//радиальная и нормальная составляющая скорости
	//double Unit_V = sqrt(nu_Earth / R_GCO);
	KA p1;



	p1.delV			= this->delV;
	p1.u			= this->u;
	p1.omeg			= this->omeg;
	p1.w			= this->w;
	p1.p			= this->p;
	p1.e			= this->e;
	p1.i			= this->i;
	p1.Time			= this->Time;
	p1.Time_start   = this->Time_start;
	p1.Time_end		= this->Time_end;

	//Переход от орбитальной СК к ГЭСК
	p1.Anom		= p1.u - p1.w;

	Vr			= sqrt(nu_Earth / p1.p)*p1.e*sin(p1.Anom);
	Vn			= sqrt(nu_Earth / p1.p)*(1 + p1.e * cos(p1.Anom));

	double R1 = (p1.p / (1 + p1.e * cos(p1.Anom)));

	p1.x		= R1 * (cos(p1.omeg)*cos(p1.u) - sin(p1.omeg)*sin(p1.u)*cos(p1.i));
	p1.y		= R1 * (sin(p1.omeg)*cos(p1.u) + cos(p1.omeg)*sin(p1.u)*cos(p1.i));
	p1.z		= R1 * sin(p1.i)*sin(p1.u);


	p1.Vx		= ((Vr * (cos(p1.u)*cos(p1.omeg) - sin(p1.u)*cos(p1.i)*sin(p1.omeg)) + Vn * (-sin(p1.u)*cos(p1.omeg) - cos(p1.i)*sin(p1.omeg)*cos(p1.u))));//(p1.delV * (Vr * (cos(p1.u)*cos(p1.omeg) - sin(p1.u)*cos(p1.i)*sin(p1.omeg)) + Vn * (-sin(p1.u)*cos(p1.omeg) - cos(p1.i)*sin(p1.omeg)*cos(p1.u))));
	p1.Vy		= ((Vr * (sin(p1.omeg) * cos(p1.u) + cos(p1.omeg) * sin(p1.u) * cos(p1.i)) + Vn * (-sin(p1.u)*sin(p1.omeg) + cos(p1.u)*cos(p1.i)*cos(p1.omeg))));//(p1.delV * (Vr * (sin(p1.omeg) * cos(p1.u) + cos(p1.omeg) * sin(p1.u) * cos(p1.i)) + Vn * (-sin(p1.u)*sin(p1.omeg) + cos(p1.u)*cos(p1.i)*cos(p1.omeg))));
	p1.Vz		= ((Vr * sin(p1.i)*sin(p1.u) + Vn * cos(p1.u)*sin(p1.i)));//(p1.delV * (Vr * sin(p1.i)*sin(p1.u) + Vn * cos(p1.u)*sin(p1.i)));

	double Mod = sqrt(pow(p1.Vx, 2) + pow(p1.Vy, 2) + pow(p1.Vz, 2));


	//Определение направля.щий косинусов
	double Cx = p1.Vx / Mod;	double Cy = p1.Vy / Mod;	double Cz = p1.Vz / Mod;

	double dVx = Cx * p1.delV;	double dVy = Cy * p1.delV;	double dVz = Cz * p1.delV;


	//Переход к БЕЗРАЗМЕРНЫМ ВЕЛИЧИНАМ
	p1.x	= p1.x / Unit_R;
	p1.y	= p1.y / Unit_R;
	p1.z	= p1.z / Unit_R;

	p1.Vx	= p1.Vx / Unit_V + dVx;
	p1.Vy	= p1.Vy / Unit_V + dVy;
	p1.Vz	= p1.Vz / Unit_V + dVz;

	Mod = sqrt(pow(p1.Vx, 2) + pow(p1.Vy, 2) + pow(p1.Vz, 2));



	p1.Time = p1.Time / Unit_T;

	//double Unit_t = Unit_R / Unit_V;

	double *xa		= new double[Time_fly * 60 * 60];
	double *ya		= new double[Time_fly * 60 * 60];
	double *za		= new double[Time_fly * 60 * 60];
	double *Vxa		= new double[Time_fly * 60 * 60];
	double *Vya		= new double[Time_fly * 60 * 60];
	double *Vza		= new double[Time_fly * 60 * 60];

	KA  p2, p3, p4, pk, k1, k2, k3, k4, kk;

	std::ofstream f;
	f.open(" report.txt");

	int n = 0; int pointer = 0;
	double h = 1.0 * 60 / Unit_T;//тут время в секундах 

	std::ofstream f4;	 f4.open(" X.txt ");	f4 << " X = [";
	std::ofstream f5;	 f5.open(" Y.txt ");	f5 << " Y = [";
	std::ofstream f6;    f6.open(" Z.txt ");	f6 << " Z = [";





	pk = p1;


	int cou = 0;
	KA e;
	double epsil = 0.0001;	double k = 8.0;

	while ((p1.Time) < (Flight_Time1 + Time_start) / Unit_T)	// Интегрирование системымы дифуров методом Рунге-Кутта 4
	{
		n++;

		//pk.print(&f);
		if (((Flight_Time1 + Time_start) / Unit_T) - pk.Time < h)
			h = ((Flight_Time1 + Time_start) / Unit_T) - pk.Time;

		xa[n] = p1.x*R_GCO;
		ya[n] = p1.y*R_GCO;
		za[n] = p1.z *R_GCO;

		Vxa[n] = p1.Vx*Unit_V;
		Vya[n] = p1.Vy*Unit_V;
		Vza[n] = p1.Vz*Unit_V;

		k1 = p1.difer(Time_end);

		p2 = p1 + k1 * (h / 2);
		k2 = p2.difer(Time_end);

		p3 = p1 + k2 * (h / 2);
		k3 = p3.difer(Time_end);

		p4 = p1 + k3 * h;
		k4 = p4.difer(Time_end);

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * (h);
		/*
		e = (k1 - k2 - k3 + k4)*(2.0 / 3.0);
		if (fabs(e.E() > epsil))
			h = h*pow(epsil / e.E(), 0.2);
		else if (fabs(e.E() < epsil / k))
			h = 2 * h;
		else if (epsil / k <= fabs(e.E()) && fabs(e.E()) <= epsil)
			h=h;*/
			//std::wcout << " h = " << h*Unit_T << std::endl;
		cou++;



		p1 = pk;
		p1.Grafic_X(&f4);	p1.Grafic_Y(&f5);	p1.Grafic_Z(&f6);

		//std::cout << pk.Time*Unit_T << std::endl;
	}

	f4 << "];";
	f5 << "];";
	f6 << "];\n plot3(X,Y,Z,'r');\n hold on;";

	double xKA = p1.x; double yKA = p1.y; double zKA = p1.z;
	int size_print = Flight_Time1 / 120.0;

	//	if (false)
	//{
	/*
	std::ofstream f2;
	f2.open(" KA_position.txt");
	f2 << " xa = [";
	for (int j = 1; j < cou; j++)
		f2 << xa[j] << ",";
	f2 << " ];\nya=[";
	for (int j = 1; j < cou; j++)
		f2 << ya[j] << ",";
	f2 << " ];\nza=[";
	for (int j = 1; j < cou; j++)
		f2 << za[j] << ",";
	f2 << "];\nplot3(xa,ya,za,'r');\nhold on;";
	f2.close();*/



	//}


	double T_L4;


	double min		= sqrt(pow((p1.x - x2 / Unit_R), 2) + pow((p1.y - y2 / Unit_R), 2) + pow((p1.z - z2 / Unit_R), 2));

	double del2Vx, del2Vy, del2Vz;		//Определение величны второго импульса скорости 

	del2Vx			= Vx2 - p1.Vx*Unit_V;	del2Vy = Vy2 - p1.Vy*Unit_V;	del2Vz = Vz2 - p1.Vz*Unit_V;

	dVx2			= del2Vx / Unit_V;	 dVy2 = del2Vy / Unit_V;	dVz2 = del2Vz / Unit_V;
	delV2			= sqrt(pow(del2Vx, 2) + pow(del2Vy, 2) + pow(del2Vz, 2));	//km/c


	std::cout << std::fixed;	std::cout.precision(5);

	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] Vxa;
	delete[] Vya;
	delete[] Vza;
	return  min;
}

void KA::Back_Rynge_Kyt()
{
	double Vr; double Vn; double R;	//радиальная и нормальная составляющая скорости
	setlocale(LC_ALL, "ru");
	KA p1;


	p1.u	= this->u;
	p1.omeg = this->omeg;
	p1.w	= this->w;
	p1.p	= this->p;
	p1.e	= this->e;
	p1.i	= this->i;
	p1.Time = this->Time;

	//p1.Time		    = this->Time_end;









	int j;
	double z_L4 = 0;
	int count = 0;
	int Start_point = 0;// Time_start / (60 * 60);
	for (int i = Start_point; i < (Month); i++)// Определение индекса массива положение точки L4
	{
		if (count < 4 + point)
		{
			if (z_L4 > zL4[i + 1] && z_L4 < zL4[i])
			{
				j = i;
				count++;
				//std::cout << " count = " << count << (H[j] + (z_L4 - zL4[j]) * (H[j + 1] - H[j]) / (zL4[j + 1] - zL4[j]))*60 << "\n\n";
				//std::cout << " Z -1 = " << zL4[i]<<" zL4 +1 = "<<zL4[i+1] << " X = " << xL4[i] << " Y = " << yL4[i]<<" h = "<<H[j]/(60.0*60.0*24) << std::endl;
			}
			else if (z_L4 == zL4[i])
				j = i;

		}
		else
			break;

	}
	count = 0;

	double T_L4;
	T_L4 = H[j] + (z_L4 - zL4[j]) * (H[j + 1] - H[j]) / (zL4[j + 1] - zL4[j]);//в минутах 


	p1.Time_end = T_L4 * 60.0;  // Конечное время должно совпадать с врменем точки L4 прохождением нисходящего узла
	p1.Time = p1.Time_end;

	double XL0; double YL0; double ZL0; //Определение координат точки L4 в момент прохождения низходящего узла орбиты
	p1.x	= xL4[j] + (T_L4 - H[j])*(xL4[j + 1] - xL4[j]) / (H[j + 1] - H[j]);
	p1.y	= yL4[j] + (T_L4 - H[j])*(yL4[j + 1] - yL4[j]) / (H[j + 1] - H[j]);
	p1.z	= zL4[j] + (T_L4 - H[j])*(zL4[j + 1] - zL4[j]) / (H[j + 1] - H[j]);

	xL		= p1.x;
	yL		= p1.y;
	zL		= p1.z;

	p1.x	= p1.x / Unit_R;
	p1.y	= p1.y / Unit_R;
	p1.z	= p1.z / Unit_R;



	double VLx = VL4X[j] + (T_L4 - H[j])*(VL4X[j + 1] - VL4X[j]) / (H[j + 1] - H[j]);
	double VLy = VL4Y[j] + (T_L4 - H[j])*(VL4Y[j + 1] - VL4Y[j]) / (H[j + 1] - H[j]);
	double VLz = VL4Z[j] + (T_L4 - H[j])*(VL4Z[j + 1] - VL4Z[j]) / (H[j + 1] - H[j]);

	VLx4 = VLx / Unit_V;
	VLy4 = VLy / Unit_V;
	VLz4 = VLz / Unit_V;




	double Mod = sqrt(pow(VLx, 2) + pow(VLy, 2) + pow(VLz, 2));
	double Cx = VLx / Mod;	double Cy = VLy / Mod;	double Cz = VLz / Mod;

	double dVx = Cx * delV3;	double dVy = Cy * delV3;	double dVz = Cz * delV3;

	p1.Vx = (VLx + dVx);		p1.Vy = (VLy + dVy);		p1.Vz = (VLz + dVz);

	p1.Vx = p1.Vx / Unit_V;
	p1.Vy = p1.Vy / Unit_V;
	p1.Vz = p1.Vz / Unit_V;

	std::cout << std::fixed;	std::cout.precision(6);
	std::cout << " ОБратное время  : Время конца =  " << p1.Time / (24.0*60.0*60.0) << " x = " << p1.x*Unit_R << " y = " << p1.y*Unit_R << " z = " << p1.z*Unit_R << " Vx = " << VLx << " Vy = " << VLy << " Vz = " << VLz << std::endl;
	p1.Time = p1.Time / Unit_T;

	//double Unit_t = Unit_R / Unit_V;

	double *xa	= new double[Time_fly * 60 * 60];
	double *ya	= new double[Time_fly * 60 * 60];
	double *za	= new double[Time_fly * 60 * 60];
	double *Vxa = new double[Time_fly * 60 * 60];
	double *Vya = new double[Time_fly * 60 * 60];
	double *Vza = new double[Time_fly * 60 * 60];

	KA  p2, p3, p4, pk, k1, k2, k3, k4, kk;

	std::ofstream f;	f.open(" report.txt");

	int n = 0; int pointer = 0;

	double h = 60 / Unit_T;//тут время в секундах 

	std::ofstream f4;	 f4.open(" X.txt ");	f4 << " X = [";
	std::ofstream f5;	 f5.open(" Y.txt ");	f5 << " Y = [";
	std::ofstream f6;    f6.open(" Z.txt ");	f6 << " Z = [";


	double R_m = 0;


	//std::cout << " ID = " << std::this_thread::get_id() << std::endl;

	pk = p1;

	int cou = 0;
	double R_max = sqrt(pow(pk.x, 2) + pow(pk.y, 2) + pow(pk.z, 2));
	double T22 = 0;
	double Z = p1.z;

	// (p1.Time >= (Flight_Time ) / Unit_T)
	do
	{
		n++;
		if (Z*Unit_R < 10.0  && Z>0)
		{
			//h = 10 / Unit_T;
		}
		else
			h = 60 / Unit_T;

		//pk.print(&f);

		xa[n]	= p1.x*R_GCO;
		ya[n]	= p1.y*R_GCO;
		za[n]	= p1.z *R_GCO;

		Vxa[n]	= p1.Vx*Unit_V;
		Vya[n]	= p1.Vy*Unit_V;
		Vza[n]	= p1.Vz*Unit_V;


		k1 = p1.Back_difer();

		p2 = p1 + k1 * (h / 2);
		k2 = p2.Back_difer();

		p3 = p1 + k2 * (h / 2);
		k3 = p3.Back_difer();

		p4 = p1 + k3 * h;
		k4 = p4.Back_difer();

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * (h);



		Z = pk.z;

		cou++;
		p1 = pk;
		//p1.Grafic_X(&f4);	p1.Grafic_Y(&f5);	p1.Grafic_Z(&f6);
		T22 += h;

	} while (Z > 0);
	p1.Time_end = T_L4 * 60.0;

	std::cout << " T2 = " << T22 * Unit_T / (24.0*60.0*60.0) << " Time_end = " << p1.Time_end / (24.0*60.0*60.0) << std::endl;
	double xKA = p1.x; double yKA = p1.y; double zKA = p1.z;
	int size_print = T2 / h;
	/*
	if (Print)
	{
		std::ofstream f2;
		f2.open(" KA_position.txt");
		f2 << " xa = [";
		for (int j = 1; j < cou; j++)
			f2 << xa[j] << ",";
		f2 << " ];\nya=[";
		for (int j = 1; j < cou; j++)
			f2 << ya[j] << ",";
		f2 << " ];\nza=[";
		for (int j = 1; j < cou; j++)
			f2 << za[j] << ",";
		f2 << "];\nplot3(xa,ya,za,'k');\nhold on;";
		f2.close();

		f2 << "];";
		f2 << "];";
		f2 << "];\n plot3(X,Y,Z,'r');\n hold on;";

		f6.close();
	}*/



	double min = 1;// sqrt(pow((XL0 / Unit_R - p1.x), 2) + pow((YL0 / Unit_R - p1.y), 2) + pow((ZL0 / Unit_R - p1.z), 2));
	//std::cout << " X = " << p1.x << " Y = " << p1.y << " Z = " << p1.z*Unit_R << std::endl;

	x = p1.x;		Vx = p1.Vx;
	y = p1.y;		Vy = p1.Vy;
	z = p1.z;		Vz = p1.Vz;

	x2 = p1.x*Unit_R;		Vx2 = p1.Vx*Unit_V;
	y2 = p1.y*Unit_R;		Vy2 = p1.Vy*Unit_V;
	z2 = p1.z*Unit_R;		Vz2 = p1.Vz*Unit_V;

	Time = p1.Time*Unit_T;

	this->T2 = T22 * Unit_T;

	Time_end = p1.Time_end;

	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] Vxa;
	delete[] Vya;
	delete[] Vza;

}

double KA::Full_Rynge_Kyt()
{
	setlocale(LC_ALL, "ru");
	double Vr; double Vn; double R;	//радиальная и нормальная составляющая скорости
	//double Unit_V = sqrt(nu_Earth / R_GCO);
	KA p1;



	p1.delV				= this->delV;
	p1.u				= this->u;
	p1.omeg				= this->omeg;
	p1.w				= this->w;
	p1.p				= this->p;
	p1.e				= this->e;
	p1.i				= this->i;
	p1.Time				= this->Time;
	p1.Time_start		= this->Time_start;


	p1.Anom = p1.u - p1.w;

	Vr = sqrt(nu_Earth / p1.p)*p1.e*sin(p1.Anom);
	Vn = sqrt(nu_Earth / p1.p)*(1 + p1.e * cos(p1.Anom));

	double R1 = (p1.p / (1 + p1.e * cos(p1.Anom)));

	p1.x = R1 * (cos(p1.omeg)*cos(p1.u) - sin(p1.omeg)*sin(p1.u)*cos(p1.i));
	p1.y = R1 * (sin(p1.omeg)*cos(p1.u) + cos(p1.omeg)*sin(p1.u)*cos(p1.i));
	p1.z = R1 * sin(p1.i)*sin(p1.u);


	p1.Vx = ((Vr * (cos(p1.u)*cos(p1.omeg) - sin(p1.u)*cos(p1.i)*sin(p1.omeg)) + Vn * (-sin(p1.u)*cos(p1.omeg) - cos(p1.i)*sin(p1.omeg)*cos(p1.u))));//(p1.delV * (Vr * (cos(p1.u)*cos(p1.omeg) - sin(p1.u)*cos(p1.i)*sin(p1.omeg)) + Vn * (-sin(p1.u)*cos(p1.omeg) - cos(p1.i)*sin(p1.omeg)*cos(p1.u))));
	p1.Vy = ((Vr * (sin(p1.omeg) * cos(p1.u) + cos(p1.omeg) * sin(p1.u) * cos(p1.i)) + Vn * (-sin(p1.u)*sin(p1.omeg) + cos(p1.u)*cos(p1.i)*cos(p1.omeg))));//(p1.delV * (Vr * (sin(p1.omeg) * cos(p1.u) + cos(p1.omeg) * sin(p1.u) * cos(p1.i)) + Vn * (-sin(p1.u)*sin(p1.omeg) + cos(p1.u)*cos(p1.i)*cos(p1.omeg))));
	p1.Vz = ((Vr * sin(p1.i)*sin(p1.u) + Vn * cos(p1.u)*sin(p1.i)));//(p1.delV * (Vr * sin(p1.i)*sin(p1.u) + Vn * cos(p1.u)*sin(p1.i)));

	double Mod = sqrt(pow(p1.Vx, 2) + pow(p1.Vy, 2) + pow(p1.Vz, 2));

	double Hh = Mod * Mod - (2 * nu_Earth) / R1;

	double Cx = p1.Vx / Mod;	double Cy = p1.Vy / Mod;	double Cz = p1.Vz / Mod;

	double dVx = Cx * p1.delV;	double dVy = Cy * p1.delV;	double dVz = Cz * p1.delV;


	//БЕЗРАЗМЕРНЫЕ ВЕЛИЧИНЫ
	p1.x		= p1.x / Unit_R;
	p1.y		= p1.y / Unit_R;
	p1.z		= p1.z / Unit_R;

	p1.Vx		= p1.Vx / Unit_V + dVx;
	p1.Vy		= p1.Vy / Unit_V + dVy;
	p1.Vz		= p1.Vz / Unit_V + dVz;

	Mod = sqrt(pow(p1.Vx, 2) + pow(p1.Vy, 2) + pow(p1.Vz, 2));
	Hh = Mod * Mod - (2 * nu_Earth) / R1;

	std::cout << std::fixed;	std::cout.precision(7);
	p1.Time = Time_start / Unit_T;

	//double Unit_t = Unit_R / Unit_V;

	double *xa = new double[Time_fly * 60 * 60];
	double *ya = new double[Time_fly * 60 * 60];
	double *za = new double[Time_fly * 60 * 60];
	double *Vxa = new double[Time_fly * 60 * 60];
	double *Vya = new double[Time_fly * 60 * 60];
	double *Vza = new double[Time_fly * 60 * 60];

	KA  p2, p3, p4, pk, k1, k2, k3, k4, kk;
	std::ofstream f;
	f.open(" report.txt");

	int n = 0; int pointer = 0;
	double h = 60.0 / Unit_T;//тут время в секундах 

	std::ofstream f4;		 f4.open(" X.txt ");		f4 << " X = [";
	std::ofstream f5;		 f5.open(" Y.txt ");		f5 << " Y = [";
	std::ofstream f6;		 f6.open(" Z.txt ");		f6 << " Z = [";
	std::ofstream pe;		 pe.open("e.txt");			pe << " e = [";
	std::ofstream pi;		 pi.open("i.txt");			pi << " i = [";
	std::ofstream pomeg;	 pomeg.open("omeg.txt");	pomeg << " omeg = [";
	std::ofstream pw;		 pw.open("w.txt");			pw << " w = [";
	std::ofstream pa;		 pa.open("a.txt");			pa << " a = [";
	std::ofstream ph;		 ph.open("h.txt");			ph << " h = [";
	std::ofstream f7;		 f7.open(" XLF.txt ");		f7 << " XL = [";
	std::ofstream f8;		 f8.open(" YLF.txt ");		f8 << " YL = [";
	std::ofstream f9;		 f9.open(" ZLF.txt ");		f9 << " ZL = [";
	std::ofstream f10;		 f10.open(" fi.txt");		f10 << " fi = [";
	std::ofstream f11;		 f11.open(" lan.txt");		f11 << " lam = [";

	double R_m = 0;

	double sx, sy, sz;

	//std::cout << " ID FULL = " << std::this_thread::get_id() << std::endl;

	pk = p1;

	p1.Time = Time_start / Unit_T;
	double R_max = sqrt(pow(pk.x, 2) + pow(pk.y, 2) + pow(pk.z, 2));
	int count = 1;


	while (p1.Time <= (Flight_Time1 + T2 + Time_start) / Unit_T)//(R_max > R_m)
	{

		if (p1.Time >= (Flight_Time1 + Time_start) / Unit_T && count == First_run)
		{
			//std::cout << "\n\n************\n*******Time = " << p1.Time*Unit_T << " TimeIm = " << Flight_Time1 + Time_start << " dvx2 = " << dVx2 << " dVy2 = " << dVy2 << " dVz2 = " << dVz2 << std::endl;

			p1.Vx = p1.Vx + dVx2;
			p1.Vy = p1.Vy + dVy2;
			p1.Vz = p1.Vz + dVz2;
			count++;
			//h = 30.0 / Unit_T;
		}
		n++;



		if (((Flight_Time1 + Time_start) / Unit_T - p1.Time) < h && ((Flight_Time1 + Time_start) / Unit_T - p1.Time) > 0)
			h = ((Flight_Time1 + Time_start) / Unit_T - p1.Time);
		else
			h = 60.0 / Unit_T;

		//if ((Flight_Time1 + T2 + Time_start) / Unit_T - p1.Time < h)
			//h = (Flight_Time1 + T2 + Time_start) / Unit_T - p1.Time;

		xa[n] = p1.x*R_GCO;
		ya[n] = p1.y*R_GCO;
		za[n] = p1.z *R_GCO;

		Vxa[n] = p1.Vx*Unit_V;
		Vya[n] = p1.Vy*Unit_V;
		Vza[n] = p1.Vz*Unit_V;


		k1 = p1.difer(Time_end);

		p2 = p1 + k1 * (h / 2);
		k2 = p2.difer(Time_end);

		p3 = p1 + k2 * (h / 2);
		k3 = p3.difer(Time_end);

		p4 = p1 + k3 * h;
		k4 = p4.difer(Time_end);

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * (h);
		/*
		int j;
		double hi = pk.Time * Unit_T;
		//std::cout << " Time = rnd = " << Time_end << "\n";
		
		for (size_t i = 0; i < Month; i++)
		{
			if (hi / (60) < H[i + 1] && hi / (60) > H[i])
			{
				j = i;
				break;
			}

			else if (hi / (60) == H[i])
			{
				j = i;
				break;

			}

		}


		//std::cout << " j = " << j << " H = " << H[j] << std::endl;


		sx = (xS[j] + (hi / (60) - H[j])*(xS[j + 1] - xS[j]) / (H[j + 1] - H[j]));
		sy = (yS[j] + (hi / (60) - H[j])*(yS[j + 1] - yS[j]) / (H[j + 1] - H[j]));
		sz = (zS[j] + (hi / (60) - H[j])*(zS[j + 1] - zS[j]) / (H[j + 1] - H[j]));

		if (Lighting(pk.x*Unit_R, pk.y*Unit_R, pk.z*Unit_R, sx, sy, sz) < 0)
			std::cout << "\n\n Dark \n";


		Lighting(pk.x*Unit_R,pk.y*Unit_R, pk.z*Unit_R, sx, sy, sz);*/


		p1 = pk;




		if (n % 9 == 0)
		{
			p1.Grafic_Orbit_Parametr(&pe, &pi, &pomeg, &pw, &pa, &ph);
			//p1.Grafic_X(&f4);	p1.Grafic_Y(&f5);	p1.Grafic_Z(&f6);
			p1.Degree(&f7, &f8, &f9);
		}

		p1.print(&f, Time_start);



	}

	f4 << "];";													f7 << "];";
	f5 << "];";													f8 << "]; \n plot(XL, YL, 'r');\n hold on;";
	f6 << "];\n plot3(X,Y,Z,'r');\n hold on;";					f10 << " ]; \n";								f11 << " ];\n plot(lam, fi,'r');";

	pe << "];";		pi << "];";		pomeg << "];";		pw << "];";		pa << "];";		ph << "];";

	double xKA = p1.x; double yKA = p1.y; double zKA = p1.z;
	int size_print = Flight_Time1 / 120.0;
	
	if (Print)
	{
		std::ofstream f2;
		f2.open(" KA_position.txt");
		f2 << " xa = [";
		for (int j = 1; j < n; j++)
			f2 << xa[j] << ",";
		f2 << " ];\nya=[";
		for (int j = 1; j < n; j++)
			f2 << ya[j] << ",";
		f2 << " ];\nza=[";
		for (int j = 1; j < n; j++)
			f2 << za[j] << ",";
		f2 << "];\nplot3(xa,ya,za,'r');\nhold on;";
		f2.close();



	}

	Mod = sqrt(pow(p1.Vx*Unit_V, 2) + pow(p1.Vy*Unit_V, 2) + pow(p1.Vz*Unit_V, 2));

	Cx = p1.Vx*Unit_V / Mod;	Cy = p1.Vy*Unit_V / Mod;	 Cz = p1.Vz*Unit_V / Mod;

	dVx = Cx * delV3;	 dVy = Cy * delV3;	 dVz = Cz * delV3;

	p1.Vx = (p1.Vx - dVx / Unit_V);		p1.Vy = (p1.Vy - dVy / Unit_V);		p1.Vz = (p1.Vz - dVz / Unit_V);


	double T_L4;


	double min = sqrt(pow((x2 / Unit_R - p1.x), 2) + pow((y2 / Unit_R - p1.y), 2) + pow((z2 / Unit_R - p1.z), 2));
	double del3Vx, del3Vy, del3Vz;
	del3Vx = Vx2 - p1.Vx*Unit_V;	del3Vy = Vy2 - p1.Vy*Unit_V;	del3Vz = Vz2 - p1.Vz*Unit_V;

	
	d = sqrt(pow(p1.x*Unit_R - xL, 2) + pow(p1.y*Unit_R - yL, 2) + pow(p1.z*Unit_R - zL, 2));
	
	double Vsum = delV + delV2 + delV3;

	std::cout << " Полное интергрирование : Время окончания =  " << p1.Time*Unit_T / (24.0*60.0*60.0) << " x = " << p1.x*Unit_R << " y = " << p1.y*Unit_R << " z = " << p1.z*Unit_R << " Vx = " << p1.Vx*Unit_V << " Vy = " << p1.Vy*Unit_V << " Vz = " << p1.Vz*Unit_V << std::endl;
	std::cout << " Xl = " << xL << " yL = " << yL << " zL = " << zL << " x = " << p1.x*Unit_R << " Y = " << p1.y*Unit_R << " z = " << p1.z*Unit_R << " Time_end = " << Time_end << " Time = " << p1.Time*Unit_T << std::endl;
	std::cout << " Промах " << d << " км " << std::endl;

	this->x = p1.x;		this->Vx = p1.Vx;
	this->y = p1.y;		this->Vy = p1.Vy;
	this->z = p1.z;		this->Vz = p1.Vz;

	this->Time = p1.Time;

	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] Vxa;
	delete[] Vya;
	delete[] Vza;
	return  Vsum;
}


void KA::Rynge_Kytta_L4()
{
	setlocale(LC_ALL, "ru");

	KA p1, p2, p3, p4, pk, k1, k2, k3, k4, kk;


	double t		= (27.332 * 24 * 60) / 6;
	double time_e	= 12 * 30 * 24 * 60.0 * 60.0;

	int c;
	for (int i = 0; i < 11*30.0*24.0; i++)
	{
		if (t<H[i + 1] && t>H[i])
			c = i;
	}
	std::cout << " H[" << c << "] = " << H[c] << " t = " << t << " H[" << c + 1 << "] = " << H[c + 1] << std::endl;


	/*
	p1.x = xL4[1] / Unit_R;	p1.Vx = VL4X[1] / Unit_V;
	p1.y = yL4[1] / Unit_R;	p1.Vy = VL4Y[1] / Unit_V;
	p1.z = zL4[1] / Unit_R;	p1.Vz = VL4Z[1] / Unit_V;*/

	p1.x = x;	p1.Vx = Vx;
	p1.y = y;	p1.Vy = Vy;
	p1.z = z;	p1.Vz = Vz;
	p1.Time = Time/Unit_T;


	double Time_dark = 0; double max_time_dark = 0;	 double max_time_dark2 = 0;	double day_dark = 0; double time_dark_midle = 0; double max_time_dark_midle = 0;	 double max_time_dark_midle_2 = 0;

	std::ofstream f4;	 f4.open(" XL4.txt ");	f4 << " XL4 = [";
	std::ofstream f5;	 f5.open(" YL4.txt ");	f5 << " YL4 = [";
	std::ofstream f6;    f6.open(" ZL4.txt ");	f6 << " ZL4 = [";
	std::ofstream f10;		 f10.open(" fi.txt");		f10 << " fi = [";
	std::ofstream f11;		 f11.open(" lan.txt");		f11 << " lam = [";
	//std::ofstream alfa;		 alfa.open(".txt");		alfa << " alfa = [";

	std::cout << " Интегрирование после всех манёвров:  Время начала =  " << p1.Time*Unit_T / (24.0*60.0*60.0) << " x = " << p1.x*Unit_R << " y = " << p1.y*Unit_R << " z = " << p1.z*Unit_R << " Vx = " << p1.Vx*Unit_V << " Vy = " << p1.Vy*Unit_V << " Vz = " << p1.Vz*Unit_V << std::endl;
	double sx, sy, sz;
	double h = 1.5 * 60.0 / Unit_T;
	int m = 0;
	while ((p1.Time) < ((time_e + Time_end) / Unit_T))
	{





		k1 = p1.difer(Time_end + time_e);

		p2 = p1 + k1 * (h / 2);
		k2 = p2.difer(Time_end + time_e);

		p3 = p1 + k2 * (h / 2);
		k3 = p3.difer(Time_end + time_e);

		p4 = p1 + k3 * h;
		k4 = p4.difer(Time_end + time_e);

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * (h);

		int j=0;
		double hi = pk.Time * Unit_T;
		//std::cout << " Time = rnd = " << Time_end << "\n";
		int End = (Time_end+time_e)/(60.0);
		for (size_t i = 0; i < Month; i++)
		{
			if (hi / (60) < H[i + 1] && hi / (60) > H[i])
			{
				j = i;
				break;
			}

			else if (hi / (60) == H[i])
			{
				j = i;
				break;

			}

		}


		//std::cout << " j = " << j << " H = " << H[j] << std::endl;
		

		sx = (xS[j] + (hi / (60) - H[j])*(xS[j + 1] - xS[j]) / (H[j + 1] - H[j])) ;
		sy = (yS[j] + (hi / (60) - H[j])*(yS[j + 1] - yS[j]) / (H[j + 1] - H[j])) ;
		sz = (zS[j] + (hi / (60) - H[j])*(zS[j + 1] - zS[j]) / (H[j + 1] - H[j])) ;

		if (Lighting(pk.x*Unit_R, pk.y*Unit_R, pk.z*Unit_R, sx, sy, sz) == light_Null)
		{
			Time_dark += h * Unit_T;
		}
		else if (Lighting(pk.x*Unit_R, pk.y*Unit_R, pk.z*Unit_R, sx, sy, sz) == light_Midle)
		{
			time_dark_midle += h * Unit_T;
		}
		/*else
		{
			Time_dark		= 0;
			time_dark_midle = 0;
		}
		max_time_dark			= Time_dark;
		max_time_dark_midle		= time_dark_midle;
		if ((max_time_dark+max_time_dark_midle) > (max_time_dark2+max_time_dark_midle_2))
		{
			max_time_dark2			= max_time_dark;
			max_time_dark_midle_2	= max_time_dark_midle;
			//std::cout << " Dark_day = " << pk.Time*Unit_T/(60.0 * 60.0 * 24.0) << std::endl;
		}*/
			


			


		p1 = pk;
		if (m % 5 == 0)
		{
			p1.Degree(&f4, &f5, &f6);

			//p1.Grafic_X(&f4);	p1.Grafic_Y(&f5);	p1.Grafic_Z(&f6);
		}
	

		//std::cout << pk.Time*Unit_T << std::endl;
	}
	std::cout << " Time = " << p1.Time / (60.0*60.0*24.0) << std::endl;;
	f4 << "];";
	f5 << "];\n plot(XL4, YL4, 'r');";
	f6 << "];\n plot3(XL4,YL4,ZL4,'r');\n hold on;";
	std::cout << " Time_Dark = " << Time_dark <<" time_dark_midle = "<< time_dark_midle <<" Full_dark_time = "<<max_time_dark2+max_time_dark_midle_2<< std::endl;
}


void KA::OptimizationSpeed()
{
	KA p1;


	p1.u = this->u;
	p1.omeg = this->omeg;
	p1.delV = this->delV;

	p1.i = this->i;				p1.Flight_Time1 = this->Flight_Time1;
	p1.p = this->p;				p1.Time = this->Time;
	p1.w = this->w;				p1.Time_end = this->Time_end;
	p1.e = this->e;				p1.Time_start = this->Time_start;



	p1.m_1		= this->m_1;
	p1.m_2		= this->m_2;
	p1.m_3		= this->m_3;
	p1.m_adap   = this->m_adap;
	p1.m_NOO	= this->m_NOO;
	p1.m_rb_D   = this->m_rb_D;

	p1.W2 = this->W2;
	p1.Y1 = this->Y1;


	p1.x2 = this->x2;			p1.Vx2 = this->Vx2;
	p1.y2 = this->y2;			p1.Vy2 = this->Vy2;
	p1.z2 = this->z2;			p1.Vz2 = this->Vz2;


	double step_u = 0.000001*M_PI / 180;
	double step_omeg = 0.000001*M_PI / 180;
	double step_delV = 0.00000001;

	//double t_V = 0.0000000000002;   double t_u = 0.0000002;	double t_omeg = 0.0000002;
	//double t_V = 0.00000000000005;   double t_u = 0.000005;	double t_omeg = 0.000005;
	double t_V = 0.00000000000001;   double t_u = 0.0000000001;	double t_omeg = 0.0000000001;
	int counter = 0;

	double df; double Mark;
	double f_step_omeg, f_step_u, f_step_delV;
	double Gradient; double norma; double f1; double f2;
	double u1, omeg1, delV1;

	double dif_U, dif_omeg, dif_delV;
	std::cout << std::fixed;		std::cout.precision(3);

	p1.metka = 0;

	double u_last, omeg_last, delV_last;



	do
	{


		f1 = p1.Rynge_Kyt();




		u_last		= p1.u;
		omeg_last	= p1.omeg;
		delV_last	= p1.delV;




		std::thread Stream1([&]() {

			p1.u  = p1.u + step_u;		f_step_u = p1.Rynge_Kyt();//std::thread Stream2([&](){f_step_u = p1.Rynge_Kyt(); });


			dif_U = (f_step_u - f1) / (step_u * 180 / M_PI);
			p1.u  = p1.u - step_u;

		});

		std::thread Stream2([&]() {

			p1.omeg  = p1.omeg + step_omeg;			f_step_omeg = p1.Rynge_Kyt();
			dif_omeg = (f_step_omeg - f1) / (step_omeg * 180 / M_PI);
			p1.omeg  = p1.omeg - step_omeg;
		});


		std::thread Stream3([&]() {
			p1.delV  = p1.delV + step_delV;			f_step_delV = p1.Rynge_Kyt();
			dif_delV = (f_step_delV - f1) / step_delV;
			p1.delV  = p1.delV - step_delV;
		});


		Stream1.join();
		Stream2.join();
		Stream3.join();



		p1.u		= p1.u - t_u * dif_U;
		p1.omeg		= p1.omeg - t_omeg * dif_omeg;
		p1.delV		= p1.delV - t_V * dif_delV;


		norma = sqrt(pow(p1.u - u_last, 2) + pow(p1.omeg - omeg_last, 2) + pow(p1.delV - delV_last, 2));

		f2 = p1.Rynge_Kyt();
		
		df = sqrt(f_step_u*f_step_u + f_step_omeg * f_step_omeg + f_step_delV * f_step_delV);
	

		if ((f2 - f1)> -t_V * 0.1*df*df)
		{
			t_V		= t_V / 1.4;
			t_omeg	= t_omeg / 1.4;
			t_u		= t_u / 1.4;
			p1.delV = p1.delV + 1.4 * t_V * dif_delV;
			p1.omeg = p1.omeg + 1.4 * t_omeg*dif_omeg;
			p1.u	= p1.u + 1.4 * t_u*dif_U;

			std::cout << " Дробление   " << " tV = " << t_V << "  t_U = " << t_u << " t_omeg = " << t_omeg << std::endl;

			counter = 0;
		}
		else
		{
			//	std::cout << " Пропускаем   " << t_V << std::endl;
			counter++;
			std::cout << " c = " << counter << " \n";
			if (counter >= 4)
			{

				t_V		= t_V * 1.1500;
				t_u		= t_u * 1.1500;
				t_omeg	= t_omeg * 1.1500;
				std::cout << " 11111" << std::endl;
				counter = 0;
			}
			//p1.u	  = p1.u    - step_metod * dif_U;
			//p1.omeg   = p1.omeg - step_metod * dif_omeg;
			//p1.delV   = p1.delV - step_metod * dif_delV;
		}


		std::cout << "\nNor = " << norma << " Res = " << fabs(f2 - f1) << " T1 = " << p1.Flight_Time1 / (24.0*60.0*60.0) << " dV1 = " << p1.delV*Unit_V << " dV2 = " << p1.delV2 << " dV3 = " << delV3 << " sdV = " << p1.delV*Unit_V + p1.delV2 + delV3 << " d = " << f2 * Unit_R << std::endl;

		//Stream1.join();
		//Stream2.join();

	} while(norma > 0.000000000001 && fabs(f2 - f1) > 0.000000000001);;
	//p1.metka = 1;
//	Mark = p1.Rynge_Kyt();
	//std::cout << " D = " << Mark * Unit_R << std::endl;
	p1.delV = p1.delV*Unit_V;
	std::cout << std::fixed;	std::cout.precision(35);
	std::cout << "\n\t\t***********Result*********** \n U = " << p1.u * 180 / M_PI << " omeg = " << p1.omeg * 180 / M_PI << " del_v = " << p1.delV << " SundV of 3 impuls = " << p1.delV + p1.delV2 + delV3 << " \n ";

	p1.delV = p1.delV / Unit_V;

	u			= p1.u;
	omeg		= p1.omeg;
	delV		= p1.delV;
	delV2		= p1.delV2;
	dVx2		= p1.dVx2;
	dVy2		= p1.dVy2;
	dVz2		= p1.dVz2;

	
}


double  KA::FunFlight(double Fly_t)
{


	Flight_Time1 = Fly_t;

	Back_Rynge_Kyt();

	Time_start	= Time - Flight_Time1;
	Time		= Time_start;

	std::cout << "Time start = " << Time_start / (60.0*60.0*24.0) << " T1 = " << Flight_Time1 / (60.0*60.0*24.0) << " Timr_end = " << Time_end << " Ful Time = " << (Time_end - Time_start) / (24.0*60.0*60.0) << std::endl;
	OptimizationSpeed();
	//Rynge_Kyt();Time

	Ful_time   = T2 + Flight_Time1;


	Time = Time_start;

	double Vs = Full_Rynge_Kyt();
	Time_end = Time_start + T2 + Flight_Time1;
	Time	 = Time_end;
	Rynge_Kytta_L4();
	/*
	std::ofstream f;	f.open("L4_position.txt");
	int start = 180 *24 - 35*24;
	int finish = 180*24;

	f << " L4x = [";
	for (int i = start; i < finish; i++)
		f << xL4[i] << ",";
	f << "];\n ";
	f << " L4y = [";

	for (int i = start; i < finish; i++)
		f << yL4[i] << ",";
	f << "];\n ";
	f << " L4z = [";
	for (int i = start; i < finish; i++)
		f << zL4[i] << ",";
	f << " ];\n plot3(L4z, L4y, L4z, 'g');\n hold on;";*/



	return   Vs;
}


double KA::MainOptimization(double dV3)
{

	double a = 18.0*24.0*60.0*60.0;
	double b = 27.0*24.0*60.0*60.0;
	//std::ofstream dd;	dd.open("dd.txt");
	delV3 = dV3;
	if (dV3 < 0.18)
	{
		a = 14.0*24.0*60.0*60.0;
		b = 19.0*24.0*60.0*60.0;
	}

	double res;
	/*
	double y, z, y0, z0;
	double fy, fz;



		y0 = a + 0.382*(b - a);
		z0 = a + b - y0;


		fy = FunFlight(y0);


		fz = FunFlight(z0);

		do
		{
			if (fy <= fz)
			{
				b = z0;
				y = a + b - y0;
				z = y0;
				fz = fy;
				fy = FunFlight(y);
			}
			else
			{
				a = y0;
				y = z0;
				z = a + b - z0;
				fy = fz;
				fz = FunFlight(z);
			}

			y0 = y;
			z0 = z;

		} while (fabs(b - a) > 0.5 * 24 * 60 * 60);
		std::cout << " Optimal Flyght time = " << ((a + b) / 2) / (24 * 60 * 60) << std::endl;
		res = ((a + b) / 2) / (24 * 60 * 60);
		OptimalT = res;*/






	res = FunFlight(25.0000746 * 24 * 60 * 60.0);//25.000946
	std::cout << " \n\n\n\n\272727%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n\n\n";
	std::ofstream f1;	f1.open(" Fly parametrs.txt");	f1 << std::fixed;	f1.precision(6);

	f1 << " T1 = " << Flight_Time1 / (24.0 * 60 * 60) << "\n T2 = " << T2 / (24.0 * 60 * 60) << "\n dV1 = " << delV * Unit_V << "\n dV2 = " << delV2 << "\n dV3 = " << delV3 << "\n sumdV = " << delV * Unit_V + delV2 + delV3 << "\n ful Time = " << (Flight_Time1 + T2) / (24.0 * 60 * 60) << "\n Start Time = " << Time_start / (24.0 * 60 * 60) << "\n End Time = " << Time_end / (24.0 * 60 * 60) << "\n u = " << u * 180.0 / M_PI << "\n omeg = " << omeg * 180.0 / M_PI << std::endl;
	return res;// FunFlight(res);

}


void Main_Calc()
{



	KA p1;



	p1.m_adap	= 120.0;// kg
	p1.m_rb_D	= 1900.0;// kg
	p1.m_NOO	= 24500.0;// kg
	p1.m_1		= 700.0;// kg
	p1.m_2		= 250.0;// kg
	p1.m_3		= 290.0; // kg
	p1.Y1		= 372.0; // c
	p1.W2		= 2989.0;	// m/c
	//p1.Time_start = 0 * 24.0*60.0*60.0;
	//p1.Time_end		= 50.0*24.0*60.0*60.0;

	p1.x = 0;		p1.Vx = 0;
	p1.y = 0;		p1.Vy = 0;
	p1.z = 0;		p1.Vz = 0;

	p1.e = 0.06811;
	p1.i = 51.0*M_PI / 180;
	p1.p = 7082.3;
	p1.w = (180.0) * M_PI / 180;
	//p1.delV3	= 0.2; //km/c




	p1.u		= 175.74692781875219793619180563837289810* M_PI / 180;
	p1.omeg		= 10.26196458688647084045442170463502407 *M_PI / 180.0;
	p1.delV		= 2.91717056233505811491113490774296224 / Unit_V;
	p1.point	= 3;

	

	p1.MainOptimization(0.2010);//0.25967
	
	//double Fly_time = 18.501699999999999999*24.0*60.0*60.0;

	//p1.MainOptimization(0.2);


	/*
	double a = 0.18;
	double b = 0.26;

	double y, z, y0, z0;
	double fy, fz;

	std::ofstream printF;	printF.open("V(N).txt");

	printF<< " V = [";

	for (int i = 0; i <= 12; i += 2)
	{

		p1.point = i;
		y0 = a + 0.382*(b - a);
		z0 = a + b - y0;


		fy = p1.MainOptimization(y0);


		fz = p1.MainOptimization(z0);
		double res;
		do
		{
			if (fy <= fz)
			{
				b = z0;
				y = a + b - y0;
				z = y0;
				fz = fy;
				fy = p1.MainOptimization(y);
			}
			else
			{
				a = y0;
				y = z0;
				z = a + b - z0;
				fy = fz;
				fz = p1.MainOptimization(z);
			}

			y0 = y;
			z0 = z;
			std::cout << " Optimal Speed V3 = " << ((a + b) / 2) << " Optimal T1 = " << p1.OptimalT << " Opt u = " << p1.u << " Opt omeg = " << p1.omeg << " Opt dv1 = " << p1.delV << " Opt dv2 = " << p1.delV2 << " Opt T2 = " << p1.T2 << " FullTime = " << p1.Ful_time << std::endl;

		} while (fabs(b - a) > 0.001);
		std::cout << "\n\n***** Optimal Speed V3 = " << ((a + b) / 2) << " Optimal T1 = " << p1.OptimalT << " Opt u = " << p1.u << " Opt omeg = " << p1.omeg << " Opt dv1 = " << p1.delV << " Opt dv2 = " << p1.delV2 << " Opt T2 = " << p1.T2 << " FullTime = " << p1.Ful_time << "****\n" << std::endl;
		res =  ((a + b) / 2);
		printF << res << ", ";
	}

	printF << " ];\n N = ";
	for (int i = 0; i <= 12; i += 2)
		printF << i << ",";
	printF << "];\n plot(N, V, 'r');\nhold on;";*/



	





	/*
	std::ofstream m;
	m.open(" L_4_position.txt");
	m << " xL = [";
	for (int i = 0; i < (Fly_time+p1.Time_start)/(60.0*60); i++)
		m << xL4[i] << ",";
	m << " ];";
	m << "\n yL = [";
	for (int i = 0; i < (Fly_time + p1.Time_start) / (60.0 * 60); i++)
		m << yL4[i] << ",";
	m << " ];";
	m << " \nzL = [";
	for (int i = 0; i < (Fly_time + p1.Time_start) / (60.0 * 60); i++)
		m << zL4[i] << ",";
	m << " ];\n plot3(xL,yL,zL, 'r');\n hold on;\n";
	m.close();*/



	//double dV = OptimizationSpeed(p1, Fly_time);
	//std::cout << " Optimal Flying time = " << Fly_time / (24.0*60.0 * 60) << " Speed impulse = " << dV << std::endl;
	/*
	std::ofstream print;	print.open(" m(T).txt");
	print << " M = [";

	double Fly_time = 3.5*24.0*60.0*60.0;//5.01699999999999999*24.0*60*60;
	do
	{
		 dV= OptimizationSpeed(p1, Fly_time);
		 print << dV << " ,";
		 Fly_time += 0.5*24.0*60.0*60.0;
	} while (Fly_time <= 6.0*24.0 * 60.0 * 60.0);
	print << "];\n Time = [";
	Fly_time = 3.5*24.0*60.0*60.0;
	do
	{

		print << Time_fly << " ,";
		Fly_time += 0.5*24.0*60.0*60.0;
	} while (Fly_time <= 6.0*24.0 * 60.0 * 60.0);
	print << "];\n plot(Time, M,'r');hold on; \n grid on;\n";*/

	//	double dV = OptimizationSpeed(p1, Fly_time);
		//std::cout << " Optimal Flying time = " << Fly_time/ (24.0*60.0*60) << " Speed impulse = " << dV << std::endl;











}




int main()
{
	Timer cloc;
	setlocale(LC_ALL, "ru");
	Main_Calc();


	if (print)
	{




	}





}











