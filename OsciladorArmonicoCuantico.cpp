// OsciladorArmonicoCuantico.cpp: define el punto de entrada de la aplicación de consola.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#define PI 3.14159265359

using namespace std;

int main(void)
{
	//declaro variables
	int k, j, N, n, nciclos;
	double lambda, sigma, x0, cte;
	double oscilaciones = 0.0; //coontador de oscilaciones
	bool paso = true; //para ver si paso por el lugar inicial

	cout << "Introduzca N: ";
	cin >> N;

	cout << "Introduzca lambda: ";
	cin >> lambda;

	//cout << "Introduzca la cte: ";
	//cin >> cte;

	x0 = N / 2.0;
	sigma = N / 16.0;

	complex <double> fi[4000], beta[4000], xi[4000], alfa[4000], b[4000], gamma[4000], k0, V[4000], s, T[4000];
	complex <double> i(0, 1), period;
	ofstream fonda("fonda.dat"), potencial("potencial.dat");

	//numero de ciclos, lambda y N
	nciclos = N / 4.0;

	//genero condiciones de contorno, alfa, s, k, V y la funcion incial
	k0 = (2 * PI * nciclos) / N;

	s = 1.0 / (4.0*(k0*k0));
	
	//relleno el vector x, primero una mitad y despues la otra
	int x[4000];
	int r = 0;
	

	for (int r = 0; r <= N; r++)
	{
		x[r] = r;
	}

	for (j = 0; j <= N; j++)
	{
		if ((x[j] > 0) || (x[j] < N))
		{
			V[j] = k0 * lambda * pow((double)(x[j]-N/2), 2.0);
			potencial << x[j] << " " << real(V[j]) << endl;
		}
		else
		{
			V[j] = 100.0;
			potencial << x[j] << " " << real(V[j]) << endl;
		}
	}
	//condiciones de contorno
	fi[0] = (0.0, 0.0);
	fi[N] = (0.0, 0.0);

	//incializo los valores de la funcion de onda
	for (j = 1; j < N; j++)
	{
		fi[j] = exp(i*k0*(double)x[j])*exp(-((x[j]-x0)*(x[j] - x0)) / (2*sigma*sigma));

		//fi[j] = exp(i*k0*(double)x[j])*exp(-8.0*((4.0*x[j] - N)*(4.0*x[j] - N)) / (N*N));
		fonda << x[j] << " " << pow(abs(fi[j]), 2.0) << " " << 0 << endl;
	}

	fonda << endl;
	fonda << endl;

	alfa[N - 1] = 0.0;

	//calculo las alfa
	for (j = N - 1; j >= 1; j--)
	{
		gamma[j] = 1.0 / (-2.0 + (2.0 * i) / s - V[j] + alfa[j]);

		alfa[j - 1] = -gamma[j];
	}


	//calculo los ciclos e inicializo las beta y las xi
	n = 1;

	beta[N - 1] = (0.0, 0.0);
	xi[N] = (0.0, 0.0);

	double modulo;
	while (n < 3000)
	{
		//pongo "pasó" en false
		paso = false;

		modulo = 0.0;

		//calculo beta con las recurrencias
		for (j = N - 1; j >= 1; j--)
		{
			beta[j - 1] = gamma[j] * (((4.0*i*fi[j]) / s - beta[j]));
		}

		//calculo los xi
		for (j = 0; j <= N; j++)
		{
			xi[j + 1] = alfa[j] * xi[j] + beta[j];
		}

		//ahora calculo fi[j+1][n]
		for (j = 0; j <= N; j++)
		{
			fi[j] = xi[j] - fi[j];
			fonda << x[j] << " " << pow(abs(fi[j]), 2.0) << " " << n << endl;
		}

		fonda << endl;
		fonda << endl;

		for (j = 1; j <= N; j++)
		{
			modulo += pow(abs(fi[j]), 2.0);
		}

		//cout << modulo << endl;

		//chequeo si pasa por el punto inicial
		/*if (x[N / 2] == x0)
		{
			oscilaciones++;

			cout << oscilaciones << endl;
		}*/
		
		//n=n+1
		n++;
	}

	//calculo el período
	period = 2.0 * PI / sqrt(2.0*k0);

	cout << period << endl;

	system("pause");

	return 0;
}

