// NaveEspacial.cpp: define el punto de entrada de la aplicación de consola.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#define PI 3.14159265359

void RungeKutta(double xi, double xi_nuevo, double h);
double f1(double Pr);
double f2(double Palfa, double r);
double f3(double Palfa, double r, double alfa, double mu, double delta, double r_, double w, double t);
double f4(double delta, double mu, double r, double r_, double alfa, double w, double t);

using namespace std;

int main(void)
{
	//declaración de constantes y reescalamiento
	double Rt, Rl, Dtl, Mt, Ml, w, m, h, G, T, V, t, H, r_, v, delta, mu, error;
	//para la nave
	double r, alfa, Pr, Palfa, x, y, vx, vy;
	//para la Luna
	double rl, xl, yl;
	double t_MAX;
	//defino variables para el Runge-Kutta de cuarto orden
	double xi[4], xi_nuevo[4]; //xi es la ec diferencial a resolver y xi_nuevo seria la nueva y calculada para el algoritmo

	//asignamos valores
	G = 6.67e-11;
	Mt = 5.9736e24;
	Ml = 0.07349e24;
	Dtl = 3.844e8;
	w = 2.6617e-6;
	Rt = 6.37816e6;
	Rl = 1.7374e6;
	m = 1000;
	delta = (G * Mt) / pow(Dtl, 3.0);
	mu = Ml / Mt;

	//asignamos la condiciones iniciales al cohete
	r = Rt / Dtl;
	v = sqrt((2.0*G*Mt)/Rt) / Dtl; //reescalada
	alfa = 0.49577543479; //cabo cañaveral
	vx = v * cos(PI / 4.0); // PI/4.0 es la posición desde la superficie a la Luna
	vy = v * sin(PI / 4.0); // PI/4.0 es la posición desde la superficie a la Luna
	Pr = v * cos(PI / 4.0 - alfa);
	Palfa = r * v * sin(PI / 4.0 - alfa);
	r_ = sqrt(1 + r * r - 2.0 * r * cos(alfa - w * 0)); //fijarse

	//ajusto el tiempo
	t = 0;
	h = 60; //en segundos, aunque reajustaremos h según se requiera
	error = pow(h, 5.0);
	t_MAX = 172800; //172800 segundos que son 2 dias
	
	//ponemos a cero las funciones (paso 1 del RK)
	xi[0] = f1(Pr);
	xi[1] = f2(Palfa, r);
	xi[2] = f3(Palfa, r, alfa, mu, delta, r_, w, t);
	xi[3] = f4(delta, mu, r, r_, alfa, w, t);

	//Evolución temporal
	while (t < t_MAX) 
	{
		RungeKutta();


		t = t + h;
	}

	return 0;
}

void RungeKutta(double xi[4], double xi_nuevo[4], double h, double Palfa, double r, double alfa, double mu, double delta, double r_, double w, double t, double Pr)
{
	double k1[4], k2[4], k3[4], k4[4];

	//k1 para f1, f2, f3 y f4
	k1[0] = h * f1(Pr);
	k1[1] = h * f2(Palfa, r);
	k1[2] = h * f3(Palfa, r, alfa, mu, delta, r_, w, t);
	k1[3] = h * f4(delta, mu, r, r_, alfa, w, t);

	//k2 para f1, f2, f3 y f4
	k2[0] = h * f1(Pr + k1[2] / 2.0);
	k2[1] = h * f2(Palfa + k1[3] / 2.0, r + k1[0] / 2.0);
	k2[2] = h * f3(Palfa + k1[3] / 2.0, r + k1[0] / 2.0, alfa + k1[1] / 2.0, mu, delta, r_, w, t+ h / 2.0);
	k2[3] = h * f4(delta, mu, r + k1[0] / 2.0, r_, alfa + k1[1] / 2.0, w, t + h / 2.0);

	//k3 para f1, f2, f3 y f4
	k3[0] = h * f1(Pr + k2[2] / 2.0);
	k3[1] = h * f2(Palfa + k2[3] / 2.0, r + k2[0] / 2.0);
	k3[2] = h * f3(Palfa + k2[3] / 2.0, r + k2[0] / 2.0, alfa + k2[1] / 2.0, mu, delta, r_, w, t + h / 2.0);
	k3[3] = h * f4(delta, mu, r + k2[0] / 2.0, r_, alfa + k2[1] / 2.0, w, t + h / 2.0);

	//k4 para f1, f2, f3 y f4
	k4[0] = h * f1(Pr + k3[2]);
	k4[1] = h * f2(Palfa + k3[3], r + k3[0]);
	k4[2] = h * f3(Palfa + k3[3], r + k3[0], alfa + k3[1], mu, delta, r_, w, t + h);
	k4[3] = h * f4(delta, mu, r + k3[0], r_, alfa + k3[1], w, t + h);

	//ahora calculo la nueva funcion
	xi_nuevo[0] = xi[0] + (1 / 6.0)*(k1[0] + 2.0*k1[1] + 2.0*k1[2] + k1[3]);
	xi_nuevo[1] = xi[1] + (1 / 6.0)*(k2[0] + 2.0*k2[1] + 2.0*k2[2] + k2[3]);
	xi_nuevo[2] = xi[2] + (1 / 6.0)*(k3[0] + 2.0*k3[1] + 2.0*k3[2] + k3[3]);
	xi_nuevo[3] = xi[3] + (1 / 6.0)*(k4[0] + 2.0*k4[1] + 2.0*k4[2] + k4[3]);

	//ahora asigno los nuevos a los viejos
	xi[0] = xi_nuevo[0];
	xi[1] = xi_nuevo[1];
	xi[2] = xi_nuevo[2];
	xi[3] = xi_nuevo[3];

	return;
}

double f1(double Pr)
{
	double dr;

	dr = Pr;

	return dr;
}

double f2(double Palfa, double r)
{
	double dalfa;

	dalfa = Palfa / (r*r);

	return dalfa;
}

double f3(double Palfa, double r, double alfa, double mu, double delta, double r_, double w, double t)
{
	double dPr;

	dPr = ((Palfa*Palfa) / pow(r, 3.0)) - delta * ((1 / (r*r) + (mu / pow(r_, 3.0))*(r - cos(alfa - w * t))));

	return dPr;
}

double f4(double delta, double mu, double r, double r_, double alfa, double w, double t)
{
	double dPalfa;

	dPalfa = -((delta*mu*r) / pow(r_, 3.0))*sin(alfa - w * t);

	return dPalfa;
}
