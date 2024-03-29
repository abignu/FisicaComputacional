// EcSchrodinger.cpp: define el punto de entrada de la aplicación de consola.
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
	double lambda;

	cout << "Introduzca N: ";
	cin >> N;

	complex <double> fi[4000],beta[4000],xi[4000],alfa[4000],b[4000],gamma[4000], k0, V[4000], s;
	complex <double> i(0, 1);
	ofstream fonda("fonda.dat"), potencial("potencial.dat");

	//numero de ciclos, lambda y N

	nciclos = N / 4.0;
	lambda = 0.3;

	//genero condiciones de contorno, alfa, s, k, V y la funcion incial
	k0 = (2*PI*nciclos) / N;

	s = 1.0/(4.0*(k0*k0));

	//adjudico los valores del potencial
	for (j = 0; j < N; j++)
	{
		if ((j < (2 * N / 5)) || (j > (3 * N / 5)))
		{
			V[j] = 0.0;
			potencial << j << " " << real(V[j]) << endl;
		}
		else
		{
			V[j] = lambda * (k0*k0);
			potencial << j << " " << real(V[j]) << endl;
		}
	}
	//condiciones de contorno
	fi[0] = (0.0, 0.0);
	fi[N] = (0.0, 0.0);

	//incializo los valores de la funcion de onda
	for (j = 1; j < N; j++)
	{
		fi[j] = exp(i*k0*(double)j)*exp(-8.0*((4.0*j-N)*(4.0*j - N))/(N*N));
		fonda << j << " " << pow(2.0, abs(fi[j])) - 1.0 << " " << 0 << endl;
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

	beta[N - 1]= (0.0, 0.0);
	xi[N] = (0.0, 0.0);
	
	double modulo;
	while (n < 3000)
	{
		modulo = 0.0;

		//calculo beta con las recurrencias
		for (j = N - 1; j >= 1; j--)
		{
			beta[j - 1]= gamma[j] * (((4.0*i*fi[j])/s - beta[j]));
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
			fonda << j << " " << pow(2.0, abs(fi[j])) - 1.0 << " " << n << endl;
		}
		fonda << endl;
		fonda << endl;

		for (j = 1; j <= N; j++)
		{
			modulo += pow(2.0, abs(fi[j]));
		}

		cout << modulo << endl;
		
		//n=n+1
		n++;
	}
	

	/*for (j = 1; j <= N; j++)
	{
		fonda << j << " " << pow(2.0, abs(fi[j])) - 1.0 << endl;
	}*/

	system("pause");

    return 0;
}

