#pragma once
void KA::print(std::ofstream *f, double Time_start)
{
	*f << std::fixed;
	f->precision(5);
	*f << "Time = " << Time * Unit_T << "( " << Time * Unit_T - Time_start << " )" << "  x = " << x * Unit_R << "\ty= " << y * Unit_R << "\tz= " << z * Unit_R << "\tVx = " << Vx * Unit_V << "\tVy = " << Vy * Unit_V << "\tVz = " << Vz * Unit_V << std::endl;
}

void KA::Grafic_X(std::ofstream *f) { *f << x * Unit_R << ", "; }
void KA::Grafic_Y(std::ofstream *f) { *f << y * Unit_R << ", "; }
void KA::Grafic_Z(std::ofstream *f) { *f << z * Unit_R << ", "; }

void KA::ShowParametr()
{
	std::cout << " X = " << x * Unit_R << "km   Y = " << y * Unit_R << "km   Z = " << z * Unit_R << "km   Vx = " << Vx * Unit_V << " km/c   Vy = " << Vy * Unit_V << " km/c   Vz = " << Vz * Unit_V << " km/c" << std::endl;
}



void KA::Grafic_Orbit_Parametr(std::ofstream *e1, std::ofstream *i1, std::ofstream *omeg1, std::ofstream *w1, std::ofstream *a1, std::ofstream *h1 )
{
	double Gx = y * Unit_R *Vz*Unit_V - z * Unit_R*Vy*Unit_V;
	double Gy = z * Unit_R*Vx*Unit_V - x * Unit_R*Vz*Unit_V;
	double Gz = x * Unit_R*Vy*Unit_V - y * Unit_R*Vx*Unit_V;

	double G = sqrt(Gx*Gx + Gy * Gy + Gz * Gz);
	double v = sqrt(pow(Vx*Unit_V, 2) + pow(Vy*Unit_V, 2) + pow(Vz*Unit_V, 2));
	double r = sqrt(pow(x*Unit_R, 2) + pow(y*Unit_R, 2) + pow(z*Unit_R, 2));

	double h = v * v - (2 * nu_Earth / r);
	double e = sqrt(1 + ((G*G)*h / (nu_Earth*nu_Earth)));

	double Gz0 = (Gz / G);
	double Gy0 = (Gy / G);
	double Gx0 = (Gx / G);

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

	double p = (G*G) / nu_Earth;
	double a = p / (1 - pow(e, 2));

	*e1 << e << ", ";
	*i1 << i * 180.0 / M_PI << ", ";
	*omeg1 << omeg * 180.0 / M_PI << ", ";
	*w1 << w * 180.0 / M_PI << ", ";
	*a1 << a << ", ";
	*h1 << Time * Unit_T / (60.0 * 60) << ", ";


	double ckolar;
	double anom, x1;
	ckolar = x * Unit_R * Vx * Unit_V + y * Unit_R * Vy * Unit_V + z * Unit_R * Vz*Unit_V;

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

	double E = 2 * atan(sqrt((1 - e) / (1 + e)*tan(anom / 2)));
	double n = sqrt(nu_Earth / pow(a, 3));
	
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
	if (MMM = First_run)
	{
		tconst = (E - e * E + n * 1500);
		Betaconst = Betta;
		lamdaconst = Betaconst * wr*tconst;
		fi = asin(z*Unit_R / r);
	}
	else
	{
		lamda = lamdaconst + (Betta - Betaconst) - wr * tconst;
		fi = asin(z*Unit_R / r);
	}
		



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

KA operator-(KA a, KA b)
{
	KA c;
	c.Time = a.Time - b.Time;

	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;

	c.Vx = a.Vx - b.Vx;
	c.Vy = a.Vy - b.Vy;
	c.Vz = a.Vz - b.Vz;
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
