#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <vector>
#include <thread>

#include "Planet_position.h"
#define ON true
#define OFF false

#include "Funcs.h"


const double RSun = 696340.0;//Sun radius
const double REarth = 6371.0;//Earth raus

#define MOON  ON
#define SUN   ON

#define First_run 1
#define Print ON


const bool print = ON;

const int Time_fly = 10 * 24; // в часах
const double step = 60.0;	// в секундах
const double eps = 100000.0;
const double R_GCO = 42164.0;
const double Time_start = 5.0 * 24.0 * 60.0 * 60.0;
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
	double Vx, Vy, Vz;
	double delV;
	double Time, d, Flight_Time;
	double V;
	int metka;
	double i, w, p, omeg, e, u;
	double Anom;
	KA difer();
	KA Start_Position();
	double Rynge_Kyt();
	void print(std::ofstream *f);
	void Grafic_X(std::ofstream *f);
	void Grafic_Y(std::ofstream *f);
	void Grafic_Z(std::ofstream *f);
	void ShowParametr();
	KA ss(KA d);

	double Y1, W2, m_adap, m_rb_D, m_ka_after_first_ipmusl, m_ka_after_second_impuls, m_fuil_for_first_impuls, m_fuil_for_second_ipmuls, m_NOO, m_1, m_2, m_3;

};

int KA::count = 0;

double Radius(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}




const double Unit_R = R_GCO;//км
const double Unit_V = sqrt(nu_Earth / Unit_R);//км/с
const double Unit_T = Unit_R / Unit_V;//c
const double Unit_A = nu_Earth / (pow(Unit_R, 2));//км/с^2
const double nu_Earth1 = 1;

KA KA::difer()
{
	count++;

	







	double Mx; double My; double Mz;
	double Sx, Sy, Sz;
	

	int j;
	double hi = Time * Unit_T;
	
		for (int i = 0; i < Time_fly; i++)
		{
			if (hi / ( 60) < H[i + 1] && hi / ( 60) > H[i])
				j = i;
			else if (hi / (60 ) == H[i])
				j = i;
		}
	
	
	//std::cout << " j = " << j << " H = " << H[j] << std::endl;
	Mx = (xM[j] + (hi / (60) - H[j])*(xM[j + 1] - xM[j]) / (H[j + 1] - H[j])) / R_GCO;
	My = (yM[j] + (hi / ( 60) - H[j])*(yM[j + 1] - yM[j]) / (H[j + 1] - H[j])) / R_GCO;
	Mz = (zM[j] + (hi / ( 60) - H[j])*(zM[j + 1] - zM[j]) / (H[j + 1] - H[j])) / R_GCO;

	Sx = ((xS[j] + (hi / (60) - H[j])*(xS[j + 1] - xS[j]) / (H[j + 1] - H[j]))) / R_GCO;
	Sy = (yS[j] + (hi / ( 60) - H[j])*(yS[j + 1] - yS[j]) / (H[j + 1] - H[j])) / R_GCO;
	Sz = (zS[j] + (hi / ( 60) - H[j])*(zS[j + 1] - zS[j]) / (H[j + 1] - H[j])) / R_GCO;



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

void KA::print(std::ofstream *f)
{
	*f << std::fixed;
	f->precision(5);
	*f << "Time = " << Time * 13713.358172732362 << "  x = " << x * R_GCO << "\ty= " << y * R_GCO << "\tz= " << z * R_GCO << "\tVx = " << Vx * 3.07466628 << "\tVy = " << Vy * 3.07466628 << "\tVz = " << Vz * 3.07466628 << std::endl;
}

void KA::Grafic_X(std::ofstream *f) { *f << x*Unit_R << ", "; }
void KA::Grafic_Y(std::ofstream *f) { *f << y*Unit_R << ", "; }
void KA::Grafic_Z(std::ofstream *f) { *f << z*Unit_R << ", "; }

void KA::ShowParametr()
{
	std::cout << " X = " << x*Unit_R << "km   Y = " << y*Unit_R << "km   Z = " << z*Unit_R << "km   Vx = " << Vx*Unit_V << " km/c   Vy = " << Vy*Unit_V << " km/c   Vz = " << Vz*Unit_V << " km/c" << std::endl;
}

KA KA::ss(KA d)
{

	return KA(d);
}

KA operator+(KA a, KA b)
{



	KA c;
	c.Time = a.Time + b.Time;

	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;

	c.Vx = a.Vx + b.Vx;
	c.Vy = a.Vy + b.Vy;
	c.Vz = a.Vz + b.Vz;


	return c;
}

KA operator*(KA a, double b)
{
	KA c;
	c.x = a.x*b;
	c.y = a.y*b;
	c.z = a.z*b;

	c.Vx = a.Vx *b;
	c.Vy = a.Vy * b;
	c.Vz = a.Vz *b;

	c.Time = a.Time * b;
	return c;
}



double Kvadrat(double x, double y, double z)
{
	return sqrt(x*x + y * y + z * z);
}


double Lighting(double x, double y, double z, double XSun, double YSun, double ZSun)
{
	double R_ka, rst, rst2, alfa, fi, beta, eta, rs, penumbra, umbra, delta, gamma;
	double Rnt = REarth;
	double Rsun = RSun;
	rs = Kvadrat(XSun, YSun, ZSun);
	R_ka = Kvadrat(x, y, z);

	alfa = acos((x*XSun + y * YSun + z * ZSun) / (R_ka*rs));
	eta = acos((Rnt + Rsun) / rs);
	delta = alfa - eta;
	umbra = 0;
	penumbra = 0;
	if (delta == M_PI / 2.0)
	{
		penumbra = 0.0;
		std::cout << " sclmdc\n";
		
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
			rst2 = Rnt / cos(delta);
			penumbra = -fabs((R_ka - rst2)*cos(delta));
		}
	}


	

	gamma = acos((Rsun - Rnt) / rs);
	fi = M_PI - gamma;
	beta = alfa - fi;
	if (beta >= 0)
	{
		rst = Rnt / acos(beta);
		umbra = (R_ka - rst)*cos(beta);
	}
	else
	{
		rst = Rnt / acos(beta);
		umbra = fabs((R_ka - rst)*cos(beta));
	}
	//std::cout << " penubmra = " << penumbra << " umbra = " << umbra << std::endl;
	if (umbra < 0)
		return umbra;
	else
		return penumbra;
}





double KA::Rynge_Kyt()
{

	double Vr; double Vn; double R;	//радиальная и нормальная составляющая скорости
	//double Unit_V = sqrt(nu_Earth / R_GCO);
	KA p1;



	p1.delV		= this->delV;
	p1.u		= this->u;
	p1.omeg		= this->omeg;
	p1.w		= this->w;
	p1.p		= this->p;
	p1.e		= this->e;
	p1.i		= this->i;
	p1.Time		= this->Time;


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
	double Cx = p1.Vx / Mod;	double Cy = p1.Vy / Mod;	double Cz = p1.Vz / Mod;

	double dVx = Cx * p1.delV;	double dVy = Cy * p1.delV;	double dVz = Cz * p1.delV;


	//БЕЗРАЗМЕРНЫЕ ВЕЛИЧИНЫ
	p1.x = p1.x / Unit_R;
	p1.y = p1.y / Unit_R;
	p1.z = p1.z / Unit_R;

	p1.Vx = p1.Vx / Unit_V + dVx;
	p1.Vy = p1.Vy / Unit_V + dVy;
	p1.Vz = p1.Vz / Unit_V + dVz;

	p1.Time = p1.Time / Unit_T;

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
	double h = 2*60 / Unit_T;//тут время в секундах 

	std::ofstream f4;	 f4.open(" X.txt ");	f4 << " X = [";
	std::ofstream f5;	 f5.open(" Y.txt ");	f5 << " Y = [";
	std::ofstream f6;    f6.open(" Z.txt ");	f6 << " Z = [";


	double R_m = 0;

	
	std::cout << " ID = " << std::this_thread::get_id() << std::endl;

	pk = p1;
	double time_dark = 0;
	double R_max = sqrt(pow(pk.x, 2) + pow(pk.y, 2) + pow(pk.z, 2));
	double sx, sy, sz;
	while (p1.Time < Flight_Time / Unit_T)//(R_max > R_m)
	{
		n++;

		//pk.print(&f);

		xa[n] = p1.x*R_GCO;
		ya[n] = p1.y*R_GCO;
		za[n] = p1.z *R_GCO;

		Vxa[n] = p1.Vx*Unit_V;
		Vya[n] = p1.Vy*Unit_V;
		Vza[n] = p1.Vz*Unit_V;
		k1 = p1.difer();

		p2 = p1 + k1 * (h / 2);
		k2 = p2.difer();

		p3 = p1 + k2 * (h / 2);
		k3 = p3.difer();

		p4 = p1 + k3 * h;
		k4 = p4.difer();

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * (h);

		if (n > First_run)
		{
			R_m = sqrt(pow(xa[n - 1] / Unit_R, 2) + pow(ya[n - 1] / Unit_R, 2) + pow(za[n - 1] / Unit_R, 2));
		}
		if (pointer == First_run)
			continue;


		int j=0;
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
		sx = (xS[j] + (hi / (60) - H[j])*(xS[j + 1] - xS[j]) / (H[j + 1] - H[j]));
		sy = (yS[j] + (hi / (60) - H[j])*(yS[j + 1] - yS[j]) / (H[j + 1] - H[j]));
		sz = (zS[j] + (hi / (60) - H[j])*(zS[j + 1] - zS[j]) / (H[j + 1] - H[j]));

		if (Lighting(pk.x*Unit_R, pk.y*Unit_R, pk.z*Unit_R, sx, sy, sz) < 0)
			time_dark += h * Unit_T;
		//R_max = sqrt(pow(pk.x, 2) + pow(pk.y, 2) + pow(pk.z, 2));

		p1 = pk;
		p1.Grafic_X(&f4);	p1.Grafic_Y(&f5);	p1.Grafic_Z(&f6);

		//std::cout << pk.Time*Unit_T << std::endl;
	}
	std::cout << " Time_dark = " << time_dark << std::endl;
	f4 << "];";
	f5 << "];";
	f6 << "];\n plot3(X,Y,Z,'r');\n hold on;";
	//std::cout << " Time_fly = " << p1.Time << std::endl;
	double xKA = p1.x; double yKA = p1.y; double zKA = p1.z;
	int size_print = Flight_Time / 120.0;
	if (Print)
	{
		std::ofstream f2;
		f2.open(" KA_position.txt");
		f2 << " xa = [";
		for (int j = 1; j < size_print; j++)
			f2 << xa[j] << ",";
		f2 << " ];\nya=[";
		for (int j = 1; j < size_print; j++)
			f2 << ya[j] << ",";
		f2 << " ];\nza=[";
		for (int j = 1; j < size_print; j++)
			f2 << za[j] << ",";
		f2 << "];\nplot3(xa,ya,za,'k');\nhold on;";
		f2.close();

		//std::ofstream f3;
		//f3.open("KAvelocity.txt");
		//for (int i = 0; i < size_print; i++)
		//	f3 << " Vx = " << Vxa[i] << " Vy = " << Vya[i] << " Vz = " << Vza[i] << std::endl;
		//f3.close();

	}


	double T_L4;


	int j;
	double z_L4 = 0;
	for (int i = 0; i < Month-100; i++)
	{
		if (z_L4 > zL4[i + 1] && z_L4 < zL4[i])
			j = i;
		else if (z_L4 == zL4[i])
			j = i;

	}

	T_L4 = H[j] + (z_L4 - zL4[j]) * (H[j + 1] - H[j]) / (zL4[j + 1] - zL4[j]);//в минутах 
	int time_L4 = T_L4;
	std::cout << " Time L4 = " << T_L4*60  << std::endl;
	double XL0; double YL0; double ZL0;
	XL0 = xL4[j] + (T_L4 - H[j])*(xL4[j + 1] - xL4[j]) / (H[j + 1] - H[j]);
	YL0 = yL4[j] + (T_L4 - H[j])*(yL4[j + 1] - yL4[j]) / (H[j + 1] - H[j]);
	ZL0 = zL4[j] + (T_L4 - H[j])*(zL4[j + 1] - zL4[j]) / (H[j + 1] - H[j]);
	//double D = sqrt(pow((XL0 - xKA), 2) + pow((YL0 - yKA), 2) + pow((ZL0 - zKA), 2));



	double VLx = VL4X[j] + (T_L4 - H[j])*(VL4X[j + 1] - VL4X[j]) / (H[j + 1] - H[j]);
	double VLy = VL4Y[j] + (T_L4 - H[j])*(VL4Y[j + 1] - VL4Y[j]) / (H[j + 1] - H[j]);
	double VLz = VL4Z[j] + (T_L4 - H[j])*(VL4Z[j + 1] - VL4Z[j]) / (H[j + 1] - H[j]);

	dVx = VLx - p1.Vx*Unit_V;	dVy = VLy - p1.Vy*Unit_V;		dVz = VLz - p1.Vz*Unit_V;
	double Modul = sqrt(pow(dVx, 2) + pow(dVy, 2) + pow(dVz, 2));
	
	V = Modul;
	//std::cout<<" |dv| = "<<Modul<<

	//std::cout << std::fixed;
	//std::cout.precision(5);
	double min = sqrt(pow((XL0 / Unit_R - p1.x), 2) + pow((YL0 / Unit_R - p1.y), 2) + pow((ZL0 / Unit_R - p1.z), 2));
	//if (p1.metka == 1)
	//{
	/*
		std::cout << " \n\nX = " << p1.x*Unit_R << " XL4 = " << XL0 << " Y = " << p1.y*Unit_R << " YL4= " << YL0 << " Z = " << p1.z*Unit_R << " ZL4 =  " << ZL0 << std::endl;
		std::cout << " Vx = " << p1.Vx*Unit_V << " VxL4 = " << VLx << " Vy = " << p1.Vy*Unit_V << " VyL4 = " << VLy << " Vz = " << p1.Vz*Unit_V << " VzL4 = " << VLz << std::endl;
		std::cout << " |Vka| = " << sqrt(pow(p1.Vx*Unit_V, 2) + pow(p1.Vy*Unit_V, 2) + pow(p1.Vz*Unit_V, 2)) << " |VL4| = " << sqrt(pow(VLx, 2) + pow(VLy, 2) + pow(VLz, 2)) << std::endl;
		std::cout << " dVx = " << dVx << " dVy = " << dVy << " dVz = " << dVz <<" |dV| = "<<Modul<< std::endl;
		std::cout << "\n\n";*/
		//}


		/*
		min = sqrt(pow((XL0 - xa[0]), 2) + pow((YL0 - ya[0]), 2) + pow((ZL0 - za[0]), 2));
		for (int j = 1; j < size_print + 1; j++)
		{
		if (sqrt(pow((XL0 - xa[j]), 2) + pow((YL0 - ya[j]), 2) + pow((ZL0 - za[j]), 2)) < min)
				min = sqrt(pow((XL0 - xa[j]), 2) + pow((YL0 - ya[j]), 2) + pow((ZL0 - za[j]), 2));
		}*/



	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] Vxa;
	delete[] Vya;
	delete[] Vza;
	return  min;
}


double OptimizationSpeed(KA p1, double Flight)
{

	p1.Flight_Time = Flight;

	double step_u = 0.000001*M_PI / 180;
	double step_omeg = 0.000001*M_PI / 180;
	double step_delV = 0.00000001;
	//double t_V = 0.000001;   double t_u = 0.01;	double t_omeg = 0.01;
	double t_V = 0.00000001;   double t_u = 0.0001;	double t_omeg = 0.0001;
	int counter = 0;

	double df; double Mark;
	double f_step_omeg, f_step_u, f_step_delV;
	double Gradient; double norma; double f1; double f2;
	double u1, omeg1, delV1;
	double X, Y, Z, Vx, Vy, Vz;
	double dif_U, dif_omeg, dif_delV;
	std::cout << std::fixed;
	std::cout.precision(3);

	p1.metka = 0;

	double u_last, omeg_last, delV_last;


	
	do
	{

	
		f1 = p1.Rynge_Kyt();




		u_last = p1.u;
		omeg_last = p1.omeg;
		delV_last = p1.delV;




		std::thread Stream1([&]() {

			p1.u = p1.u + step_u;		f_step_u = p1.Rynge_Kyt();//std::thread Stream2([&](){f_step_u = p1.Rynge_Kyt(); });


			dif_U = (f_step_u - f1) / (step_u * 180 / M_PI);
			p1.u = p1.u - step_u;

		});

		std::thread Stream2([&]() {

			p1.omeg = p1.omeg + step_omeg;			f_step_omeg = p1.Rynge_Kyt();
			dif_omeg = (f_step_omeg - f1) / (step_omeg * 180 / M_PI);
			p1.omeg = p1.omeg - step_omeg;
		});


		std::thread Stream3([&]() {
			p1.delV = p1.delV + step_delV;			f_step_delV = p1.Rynge_Kyt();
			dif_delV = (f_step_delV - f1) / step_delV;
			p1.delV = p1.delV - step_delV;
		});


		Stream1.join();
		Stream2.join();
		Stream3.join();

		//std::cout << "dif_U = " << dif_U << " dif_omeg = " << dif_omeg << " dif_V = " << dif_delV << " norma = " << norma << std::endl;

		p1.u = p1.u - t_u * dif_U;
		p1.omeg = p1.omeg - t_omeg * dif_omeg;
		p1.delV = p1.delV - t_V * dif_delV;


		norma = sqrt(pow(p1.u - u_last, 2) + pow(p1.omeg - omeg_last, 2) + pow(p1.delV - delV_last, 2));//sqrt(pow(dif_U, 2) + pow(dif_omeg, 2)+pow(dif_delV,2));// +pow(dif_delV, 2));

		f2 = p1.Rynge_Kyt();
		//X = p1.x;	Y = p1.y;	Z = p1.z;	Vx = p1.Vx;	 Vy = p1.Vy;	Vz = p1.Vz;
		df = sqrt(f_step_u*f_step_u + f_step_omeg * f_step_omeg + f_step_delV * f_step_delV);
		std::cout << " U = " << p1.u * 180 / M_PI << " omeg = " << p1.omeg * 180 / M_PI << " del_v = " << p1.delV *Unit_V<<" dv2 = "<<p1.V<< " d = " << f2 * Unit_R << std::endl;

		if ((f2 - f1) > -t_V * 0.1*df*df)
		{
			t_V = t_V / 1.8;
			t_omeg = t_omeg / 1.8;
			t_u = t_u / 1.8;
			p1.delV = p1.delV + 2 * t_V * dif_delV;
			p1.omeg = p1.omeg + 2 * t_omeg*dif_omeg;
			p1.u = p1.u + 2 * t_u*dif_U;
			std::cout << " Дробление   " << " tV = " << t_V << "  t_U = " << t_u << " t_omeg = " << t_omeg << std::endl;

			counter = 0;
		}
		else
		{
			std::cout << " Пропускаем   " << t_V << std::endl;
			counter++;
			if (counter > 8)
			{

				t_V = t_V * 1.15;
				t_u = t_u * 1.15;
				t_omeg = t_omeg * 1.15;
			}
			//p1.u	  = p1.u    - step_metod * dif_U;
			//p1.omeg   = p1.omeg - step_metod * dif_omeg;
			//p1.delV   = p1.delV - step_metod * dif_delV;
		}


		std::cout << " Norma = " << norma << " Res = " << fabs(f2 - f1) << " dv" << p1.V << " sumdV = " << p1.delV*Unit_V + p1.V << std::endl;

		//Stream1.join();
		//Stream2.join();

	} while (norma > 0.0000001 || fabs(f2 - f1) > 0.0000001);
	p1.metka = 1;
	Mark = p1.Rynge_Kyt();
	std::cout << " D = " << Mark * Unit_R << std::endl;
	std::cout << "\n\t\t***********Result*********** \n U = " << p1.u * 180 / M_PI << " omeg = " << p1.omeg * 180 / M_PI << " del_v = " << p1.delV*Unit_V << " \n ";


	
	p1.delV = p1.delV*Unit_V;

	p1.m_NOO = p1.m_NOO - p1.m_1 - p1.m_2 - p1.m_3;

	p1.m_ka_after_first_ipmusl = p1.m_NOO*exp(-p1.delV*1000 / (p1.Y1*g));

	p1.m_fuil_for_first_impuls = p1.m_NOO - p1.m_ka_after_first_ipmusl;

	p1.m_ka_after_first_ipmusl = p1.m_ka_after_first_ipmusl - p1.m_adap - p1.m_rb_D;


	p1.m_ka_after_second_impuls = p1.m_ka_after_first_ipmusl*exp(-p1.V * 1000 / p1.W2);
	p1.m_fuil_for_second_ipmuls = p1.m_ka_after_first_ipmusl - p1.m_ka_after_second_impuls;

	std::cout <<"\n\n***************************************\n Time = "<<p1.Flight_Time<< " Massa KA = " << p1.m_ka_after_second_impuls<<" M fuil for 1-st impuls = "<<p1.m_fuil_for_first_impuls<<" m fuil for 2-st impuls = "<<p1.m_fuil_for_second_ipmuls<<"\n***********************************8\n\n"<< std::endl << std::endl;




	return p1.m_ka_after_second_impuls;// p1.V + p1.delV;
}

void Main_Calc()
{







	KA p1;



	p1.m_adap = 120.0;// kg
	p1.m_rb_D = 1900.0;// kg
	p1.m_NOO = 24500.0;// kg
	p1.m_1 = 700.0;// kg
	p1.m_2 = 250.0;// kg
	p1.m_3 = 290.0; // kg
	p1.Y1 = 372.0; // c
	p1.W2 = 2989.0;	// m/c


	p1.x = 0;		p1.Vx = 0;
	p1.y = 0;		p1.Vy = 0;
	p1.z = 0;		p1.Vz = 0;

	p1.e = 0.000006811;
	p1.i = 51.0*M_PI / 180;
	p1.p = 7082.3;
	p1.w = (0.0) * M_PI / 180;

	p1.Time = 0;

	

	p1.u = 355.9807466069880774739 * M_PI / 180;
	p1.omeg = 12.8790132656588021121 * M_PI / 180;
	p1.delV = 3.0030258268406657685 /Unit_V; 

	double Fly_time =4.01* 24.0*60.0*60.0;
	p1.Flight_Time = Fly_time;
	//p1.Rynge_Kyt();

	//double dV = OptimizationSpeed(p1, Fly_time);
	//std::cout << " Optimal Flying time = " << Fly_time / (24.0*60.0 * 60) << " Speed impulse = " << dV << std::endl;

	p1.Rynge_Kyt();
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



	/*
	double a = 4.5*24.0*60.0*60.0;
	double b = 5.5*24.0*60.0*60.0;

	double y, z, y0, z0;
	double fy, fz;

	y0 = a + 0.382*(b - a);
	z0 = a + b - y0;

	
	fy = OptimizationSpeed(p1, y0);

	
	fz = OptimizationSpeed(p1, z0);

	do
	{
		if (fy <= fz)
		{
			b = z0;
			y = a + b - y0;
			z = y0;
			fz = fy;
			fy = OptimizationSpeed(p1, y);
		}
		else
		{
			a = y0;
			y = z0;
			z = a + b - z0;
			fy = fz;
			fz = OptimizationSpeed(p1, z);
		}

		y0 = y;
		z0 = z;

	} while (fabs(b - a) > 0.1 * 24 * 60 * 60);
	std::cout << " Optimal Flyght time = " << ((a + b) / 2)/(24*60*60) << std::endl;*/
	
	
	


	
	
}

















int main()
{
	Timer cloc;
	setlocale(LC_ALL, "ru");
	Main_Calc();

	
	
		std::ofstream m;
		m.open(" L4_position.txt");
		m << " xL = [";
		for (int i = 0; i < Time_fly+18*24; i++)
			m << xL4[i] << ",";
		m << " ];";
		m << "\n yL = [";
		for (int i = 0; i < Time_fly + 18 * 24; i++)
			m << yL4[i] << ",";
		m << " ];";
		m << " \nzL = [";
		for (int i = 0; i < Time_fly + 18 * 24; i++)
			m << zL4[i] << ",";
		m << " ];\n plot3(xL,yL,zL, 'r');\n hold on;\n";
		m.close();

	
	





}

