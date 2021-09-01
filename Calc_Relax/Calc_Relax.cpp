

#include <iostream>
#include <fstream>
#include "Header1.h"



class KA
{
public:

	double wx, wy, wz;
	double y, fi, theta;
	double Time;
	double omeg_X, omeg_Y, omeg_Z;
	KA difer();
	void printPerametrs(std::ofstream *f);
	
};

KA operator+(KA a, KA b)
{

	KA c;

	c.Time = a.Time + b.Time;

	c.wx = a.wx + b.wx;
	c.wy = a.wy + b.wy;
	c.wz = a.wz + b.wz;
	
	c.y			= a.y	  +		b.y;
	c.fi		= a.fi    +		b.fi;
	c.theta		= a.theta +		b.theta;

	c.omeg_X = a.omeg_X + b.omeg_X;
	c.omeg_Y = a.omeg_Y + b.omeg_Y;
	c.omeg_Z = a.omeg_Z + b.omeg_Z;

	return c;
}

KA operator-(KA a, KA b)
{
	KA c;
	c.Time = a.Time - b.Time;

	c.wx = a.wx - b.wx;
	c.wy = a.wy - b.wy;
	c.wz = a.wz - b.wz;

	c.y			= a.y		- b.y;
	c.fi		= a.fi		- b.fi;
	c.theta     = a.theta	- b.theta;

	c.omeg_X = a.omeg_X - b.omeg_X;
	c.omeg_Y = a.omeg_Y - b.omeg_Y;
	c.omeg_Z = a.omeg_Z - b.omeg_Z;

	return c;
}

KA operator*(KA a, double b)
{
	KA c;
	c.wx = a.wx*b;
	c.wy = a.wy*b;
	c.wz = a.wz*b;

	c.y		= a.y*b;
	c.fi    = a.fi*b;
	c.theta = a.theta*b;

	c.omeg_X = a.omeg_X*b;
	c.omeg_Y = a.omeg_Y*b;
	c.omeg_Z = a.omeg_Z*b;

	c.Time = a.Time * b;
	return c;
}

void KA::printPerametrs(std::ofstream *f)
{
	if(!Radian)
		*f << " Time = " << Time << " wx = " << wx  << " wy = " << wy  << " wz = " << wz  << " y = " << y  << " fi = " << fi  << " theta = " << theta<<" omegX = "<<omeg_X<<" omegY = "<<omeg_Y <<" omegZ = "<<omeg_Z<< std::endl;
	else
		*f << " Time = " << Time << " wx = " << wx * 180 / myData::M_PI << " wy = " << wy * 180 / myData::M_PI << " wz = " << wz * 180 / myData::M_PI << " y = " << y * 180 / myData::M_PI << " fi = " << fi * 180 / myData::M_PI << " theta = " << theta * 180 / myData::M_PI << (wx-omeg_X) * 30 / myData::M_PI << " omegY = " << (wy-omeg_Y) * 30 / myData::M_PI << " omegZ = " << (wz-omeg_Z) * 30 / myData::M_PI << std::endl;

}

double ControlMomentX(double y, double wx)     { 
	double k = -myData::Kpx*y - myData::Krx*wx;
	return   k;
}

double ControlMomentY(double fi, double wy)     { 
	double k = -myData::Kpy*fi - myData::Kry*wy;
	return  k;
}

double ControlMomentZ( double theta, double wz)  { 
	double k = -myData::Kpy*theta - myData::Krx*wz;
	return k;
}


KA KA::difer()
{
	KA d;



	d.wx		= (1 / myData::J1)*(ControlMomentX(y, wx) + wz *myData::Jm*omeg_Y  - wy*myData::Jm*omeg_Z  - (myData::J3 - myData::J2)*wz*wy);    
	d.wy		= (1 / myData::J2)*(ControlMomentY(fi, wy) +  wx* myData::Jm*omeg_Z - wz * myData::Jm*omeg_X - (myData::J1 - myData::J3)*wx*wz);
	d.wz		= (1 / myData::J3)*(ControlMomentZ(theta, wz) + wy * myData::Jm*omeg_X - wx * myData::Jm*omeg_Y - (myData::J2 - myData::J1)*wx*wy);

	d.y			= wx;
	d.fi		= wy;
	d.theta		= wz;


	d.omeg_X	= -ControlMomentX(y, wx)/myData::Jm;
	d.omeg_Y	= -ControlMomentY(fi, wy) / myData::Jm;
	d.omeg_Z	= -ControlMomentZ(theta, wz) / myData::Jm;

	d.Time		= 1;

	return d;

}

int main()
{
	
	KA p1, p2, p3, p4, pk, k1, k2, k3, k4, kk;
	double h = 0.80;

	p1.Time = 0;

	p1.y		= myData::y_start;
	p1.fi		= myData::fi_start;
	p1.theta	= myData::theta_start;

	p1.wx		= myData::wx_start;
	p1.wy		= myData::wy_start;
	p1.wz		= myData::wz_start;

	p1.omeg_X	= myData::wx_start;
	p1.omeg_Y	= myData::wy_start;
	p1.omeg_Z	= myData::wz_start;

	int size = myData::Time_end / h;

	double *Wx		= new double[size+1];
	double *Wy		= new double[size+1];
	double *Wz		= new double[size+1];
	double *Y		= new double[size+1];
	double *Fi		= new double[size+1];
	double *Theta	= new double[size+1];
	double *Ox		= new double[size+1];
	double *Oy		= new double[size+1];
	double *Oz		= new double[size+1];

	std::ofstream f1; f1.open(" Param.txt");
	std::ofstream f2;	f2.open(" Matlab W.txt");
	std::ofstream f3;	f3.open(" Matlab Deg.txt");
	std::ofstream f4;	f4.open(" Matlab Omeg.txt");
	std::ofstream f5;	f5.open(" Matlab M.txt");
	
	int i = 0;
	pk = p1;

	double max_X = 0, max_Y = 0, max_Z = 0;
	double max_Mx = 0, max_My = 0, max_Mz = 0;
	int metka = 1;
	while (p1.Time < myData::Time_end)
	{
		

		//pk.print(&f);
		//if ((myData::Time_end - pk.Time) < h)
			//h = (myData::Time_end - pk.Time) - pk.Time;

		p1.printPerametrs(&f1);

		k1 = p1.difer();

		

		p2 = p1 + k1 * (h / 2);

		

		k2 = p2.difer();

		p3 = p1 + k2 * (h / 2);
		k3 = p3.difer();

		p4 = p1 + k3 * h;
		k4 = p4.difer();

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * (h);
	

		Wx[i]		= pk.wx;
		Wy[i]		= pk.wy;
		Wz[i]		= pk.wz;
		Y[i]		= pk.y;
		Fi[i]		= pk.fi;
		Theta[i]	= pk.theta;
		Ox[i]		= pk.omeg_X ;
		Oy[i]		= pk.omeg_Y ;
		Oz[i]		= pk.omeg_Z;
		//std::cout << " Max  OmegX= " << (-Wx[i] + Ox[i]) * 30 / myData::M_PI << std::endl;


		if (max_X < (-Wx[i] + Ox[i]) * 30 / myData::M_PI)
			max_X = (-Wx[i] + Ox[i]) * 30 / myData::M_PI;

		if (max_Y < (-Wy[i] + Oy[i]) * 30 / myData::M_PI)
			max_Y = (-Wy[i] + Oy[i]) * 30 / myData::M_PI;

		
		if (max_Z < (-Wz[i] + Oz[i]) * 30 / myData::M_PI)
			max_Z = (-Wz[i] + Oz[i]) * 30 / myData::M_PI;
		

		if (max_Mx < fabs(ControlMomentX(Y[i], Wx[i])))
			max_Mx = fabs(ControlMomentX(Y[i] , Wx[i]) );

		if (max_My < fabs(ControlMomentY(Fi[i], Wy[i])))
			max_My = fabs(ControlMomentY(Fi[i] , Wy[i]));



		if (max_Mz < fabs(ControlMomentZ(Theta[i], Wz[i] )))
			max_Mz = fabs(ControlMomentZ(Theta[i], Wz[i]));

			i++;
		p1 = pk;


		if ((pk.wx * 180 / myData::M_PI < myData::wx_min) && (pk.wy * 180.0 / myData::M_PI < myData::wy_min) && (pk.wz * 180.0 / myData::M_PI < myData::wz_min) && metka==First_Run)
		{
			std::cout << " Time relax = " << p1.Time << " Wx = " << pk.wx * 180 / myData::M_PI << " Wy = " << pk.wy * 180 / myData::M_PI << " Wz = " << pk.wz * 180 / myData::M_PI << std::endl;
			metka++;
		}
			

		//p1.y		= pk.wx;
		//p1.fi		= pk.wy;
		//p1.theta	= pk.wz;
		
	}
	//max_X = (-Wx[0] + Ox[0]) * 30 / myData::M_PI;
	//max_Y = (-Wy[0] + Oy[0]) * 30 / myData::M_PI;
	//max_Z = (-Wz[0] + Oz[0]) * 30 / myData::M_PI;

	for (int i = 1; i < size+100; i++)
	{
		

	}
	std::cout << "\n Max_wX = " << max_X << " Max_wY = " << max_Y << " Max_wZ = " << max_Z << " Max_N = " << myData::N_max << std::endl;

	std::cout << "\n Max_Mx = " << max_Mx*1000 << " Max_My = " << max_My*1000 << " Max_Mz = " << max_Mz*1000 << " Max_M = " << myData::M_max << std::endl;
	

	f2 << " wx =[";
	for (int i = 0; i < size; i++)
		f2 << Wx[i] * 180 / myData::M_PI << ",";
	f2 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f2 << i * h << ",";
	f2 << "];\n plot(Time, wx,'r');\n hold on;\n";


	f2 << "\n wy =[";
	for (int i = 0; i < size; i++)
		f2 << Wy[i] * 180 / myData::M_PI << ",";
	f2 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f2 << i * h << ",";
	f2 << "];\n plot(Time, wy,'r');\n hold on;\n";

	f2 << " \nwz =[";
	for (int i = 0; i < size; i++)
		f2 << Wz[i] * 180 / myData::M_PI << ",";
	f2 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f2 << i * h << ",";
	f2 << "];\n plot(Time, wz,'r');\n hold on;\n";


	f3 << "\n y =[";
	for (int i = 0; i < size; i++)
		f3 << Y[i] * 180 / myData::M_PI << ",";
	f3 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f3 << i * h << ",";
	f3 << "];\n plot(Time, y,'r');\n hold on;\n";


	f3 << "\n fi =[";
	for (int i = 0; i < size; i++)
		f3 << Fi[i] * 180 / myData::M_PI << ",";
	f3 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f3 << i * h << ",";
	f3 << "];\n plot(Time, fi,'r');\n hold on;\n";



	f3 << "\n theta =[";
	for (int i = 0; i < size; i++)
		f3 << Theta[i] * 180 / myData::M_PI << ",";
	f3 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f3 << i * h << ",";
	f3 << "];\n plot(Time, theta,'r');\n hold on;\n";


	f4 << "\n omegx =[";
	for (int i = 0; i < size; i++)
	{
		f4 << (-Wx[i] + Ox[i]) * 30 / myData::M_PI << ",";
		//if((Wx[i] - Ox[i]))
	}
		
	f4 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f4 << i * h << ",";
	f4 << "];\n plot(Time, omegx,'r');\n hold on;\n";


	f4 << "\n omegy =[";
	for (int i = 0; i < size; i++)
		f4 << (-Wy[i] + Oy[i]) * 30 / myData::M_PI << ",";
	f4 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f4 << i * h << ",";
	f4 << "];\n plot(Time, omegy,'r');\n hold on;\n";


	f4 << "\n omegz =[";
	for (int i = 0; i < size; i++)
		f4 << (-Wz[i] + Oz[i]) * 30 / myData::M_PI << ",";
	f4 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f4 << i * h << ",";
	f4 << "];\n plot(Time, omegz,'r');\n hold on;\n";


	f5 << "\n Mx =[";
	for (int i = 0; i < size; i++)
		f5 << ControlMomentX(Y[i], Wx[i])*1000 << ",";
	f5 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f5 << i * h << ",";
	f5 << "];\n plot(Time, Mx,'r');\n hold on;\n";


	f5 << "\n My =[";
	for (int i = 0; i < size; i++)
		f5 << ControlMomentY(Fi[i], Wy[i])*1000 << ",";
	f5 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f5 << i * h << ",";
	f5 << "];\n plot(Time, My,'r');\n hold on;\n";


	f5 << "\n Mz =[";
	for (int i = 0; i < size; i++)
		f5 << ControlMomentZ(Theta[i], Wz[i]) *1000<< ",";
	f5 << "];\n Time = [";

	for (int i = 0; i < size; i++)
		f5 << i * h << ",";
	f5 << "];\n plot(Time, Mz,'r');\n hold on;\n";

	delete[] Wx;
	delete[] Wy;
	delete[] Wz;
	delete[] Y;
	delete[] Fi;
	delete[] Theta;
	delete[] Ox;
	delete[] Oy;
	delete[] Oz;

}

