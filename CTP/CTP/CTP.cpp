



#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <vector>
#include <thread>
const bool print = 0;

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
void Koordinats(double *x, double *y, double *z, double *h, int size)
{
	Timer a;
	double step;//в часах
	std::ifstream fin;
	int count = 0;
	std::string name = "horizons.txt";
	fin.open(name);
	std::string line;
	if (!fin.is_open())
		std::cout << " Mistake " << std::endl;
	else
		std::cout << "Fail is open " << std::endl;

	std::string lineMark;
	std::string mark = "";



	int j1 = 0;
	while (std::getline(fin, line))//!fin.eof()
	{
		std::istringstream iss(line);
		iss >> lineMark;
		if (lineMark == "Step-size")
			iss >> mark >> step;
		if (lineMark == "$$SOE")
		{
			while (mark != "$$EOE")
			{
				fin >> mark;
				if (mark == "TDB")
				{
					if (count < size)
					{
						fin >> x[count] >> y[count] >> z[count];
						count++;
					}
					else
						continue;
				}

			}
		}
	}
	for (int i = 0; i < size; i++)
		h[i] = step * i;
	fin.close();
}
void Text()
{

	int size = 24 * 60 * 30;
	double *x = new double[size];
	double *y = new double[size];
	double *z = new double[size];
	double *h = new double[size];


	Koordinats(x, y, z, h, size);

	double hi = 15.865 * 60;



	double  xi, yi, zi;
	int j;
	for (int i = 0; i < size; i++)
	{
		if (hi<h[i + 1] && hi>h[i])
			j = i;
		else if (hi == h[i])
			j = i;

	}

	xi = x[j] + (hi - h[j])*(x[j + 1] - x[j]) / (h[j + 1] - h[j]);
	yi = y[j] + (hi - h[j])*(y[j + 1] - y[j]) / (h[j + 1] - h[j]);
	zi = z[j] + (hi - h[j])*(z[j + 1] - z[j]) / (h[j + 1] - h[j]);
	std::cout << x[j] << " " << " |h = " << hi << " x= " << xi << "| " << x[j + 1] << " \t";
	std::cout << y[j] << " " << " |h = " << hi << " y= " << yi << "| " << y[j + 1] << " \t";
	std::cout << z[j] << " " << " |h = " << hi << " z= " << zi << "| " << z[j + 1] << " \n";

	if (print)
	{
		std::cout << " x= ";
		for (int i = 0; i < size; i++)
			std::cout << x[i] << " ";
		std::cout << std::endl;
		std::cout << " y= ";
		for (int i = 0; i < size; i++)
			std::cout << y[i] << " ";
		std::cout << std::endl;
		std::cout << " z= ";
		for (int i = 0; i < size; i++)
			std::cout << z[i] << " ";
		std::cout << "\n\n\n";

	}


	delete[] x;
	delete[] y;
	delete[] z;
	delete[] h;

}

const double M_PI = 3.1415926535;
const double nu = 398600.0;
const double ed = 0.85;
const double eevti = 0.05;
const double cc = 880.0;
const double mas = 2000.0;
const double H_orb = 378629.0;
const double r0 = 6371.0;
const double alfa = 5.0*M_PI / 180.0;
const double Fs = 5.51;
const double umbra = 286.0*60.0;
const double beta = M_PI - 0.5*umbra*sqrt((nu / (pow(r0 + H_orb, 3))));
const double G = 5.68 / 100000000.0;

class Test
{
public:
	double  x, y, z;
	Test difer();
	void print();
	double E();
	void K3(Test k1, Test k2);
	void K4(Test k1, Test k2, Test k3);
	void K5(Test k1, Test k2, Test k3, Test k4);
	void k6(Test k1, Test k2, Test k3, Test k4, Test k5);

};

double c1(double gama)
{
	return (cos(gama)*cos(alfa) + fabs(cos(gama)*cos(alfa))) / 2.0;
}
double c2(double gama)
{
	return (-cos(gama)*cos(alfa) + fabs(-cos(gama)*cos(alfa))) / 2.0;
}
double c3(double gama)
{
	return (-sin(gama)*cos(alfa) + fabs(-sin(gama)*cos(alfa))) / 2.0;
}
double c4(double gama)
{
	return (sin(gama)*cos(alfa) + fabs(sin(gama)*cos(alfa))) / 2.0;
}
double c5(double gama)
{
	return (-sin(alfa) + fabs(-sin(alfa))) / 2.0;
}
double c6(double gama)
{
	return (sin((alfa)+fabs(sin(alfa)))) / 2.0;
}

double Qs(double gamma)
{
	if (gamma > beta && gamma < (2 * M_PI - beta))
		return 0;
	else
		return 1367.6;
}


Test Test::difer()
{
	Test d;

	//d.z = (2 * h*z) / (h*h + 1);
	d.y = ((((ed*Fs*(c5(x)*c6(x))) + (eevti * Fs*(c1(x) + c2(x) + c3(x) + c4(x)))))*Qs(x) + 4255.0) / (cc*mas) - (Fs*G * 2 * ed*pow(y, 4)) / (cc*mas);
	d.x = 1;
	return d;
}

void Test::print()
{
	std::cout << " x = " << x << " y= " << y-273 << " z= " << z << std::endl;
}

double Test::E()
{
	if (z > y)
		return z;
	else
		return y;
}

Test operator+(Test a, Test b)
{
	Test c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return c;
}


Test operator-(Test a, Test b)
{
	Test c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}


Test operator*(Test a, double b)
{
	Test c;
	c.x = a.x*b;
	c.y = a.y*b;
	c.z = a.z*b;
	return c;
}

double funx(double x)
{
	return 4 * exp(-x) - exp(2 * x);
}
double funy(double x)
{
	return   exp(2 * x);
}
void Fun()
{
	std::cout << " Fun \n\n";
	double x = 0, step = 0.1, y = 0;
	while (x <= 1)
	{
		y = x * sin(x) + cos(x);
		std::cout << " yFun = " << y << " x = " << x << std::endl;
		x += step;
	}
}





void Rynge_Kyt()
{
	using namespace std;
	Test p1, p2, p3, p4, pk, k1, k2, k2k, k3, k3k, k4, kk, p5, k5;

	p1.x = 0;
	p1.y = 293.0;
	p1.z = 3;
	double h = 1;
	double x = 0;
	
	
	double right_board = 2377400.312;
	double eps = 0.01;
	int count = 0;
	do
	{
		p1.print();
		
		k1 = p1.difer();

		p2 = p1 + k1 * (h / 2.0);
		k2 = p2.difer();

		p3 = p1 + k2 * (h / 2.0);
		k3 = p3.difer();

		p4 = p1 + k3 * h;
		k4 = p4.difer();

		kk = (k1 + k2 * 2 + k3 * 2 + k4)*(1.0 / 6.0);
		pk = p1 + kk * h;
		p1 = pk;

		

		x += h;

		

		if (right_board - x < h)
			h = right_board - x;


		count++;


	} while (x != right_board);
	std::cout << "count = " << count << std::endl;
}





int main()
{
	Rynge_Kyt();
	//Fun();
	//Text();
	double  gamma = 0;
	std::ofstream f;	f.open("Q sun.txt");
	f << " Q = [";
	double Q_sun = 0;
	do
	{
		Q_sun = (((ed*Fs*(c5(gamma)*c6(gamma))) + (eevti * Fs*(c1(gamma) + c2(gamma) + c3(gamma) + c4(gamma)))))*Qs(gamma)+4255.0;
		f << Q_sun << ", ";
		std::cout << " gamma = " << gamma * 180.0 / M_PI << " Qs = " << Q_sun << std::endl;
		gamma += 0.2*M_PI / 180.0;
	
	} while (gamma < 360.0*M_PI / 180.0);

	f << " ];\n gamma = [";
	gamma = 0;
	do
	{
			f << gamma*180.0/M_PI << ", ";

		gamma += 0.2*M_PI / 180.0;
	} while (gamma < 360.0*M_PI / 180.0);

	f << "];\n plot(gamma, Q,'r');\n hold on;\n grid on;";
}


