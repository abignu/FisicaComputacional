#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <cstdlib>

//#include "aleatorios2.h"

void MuestraMatriz(int s[100][100], int N);
double Energia(int s[100][100], int N);
double minimo(double a, double b);
double magnetizacion(int N, int s[100][100]);

using namespace std;

int main(void)
{
    //defino variables
    int N,i,j,n,m; // N es la dimension de la matriz, n y m son las coordenadas elegidas al azar
    int s[100][100];
    double P,T,hi,beta,E,dE; // hi es el numero aleatoriom P es la probabilidad, T la temperatura y E es la energia
    double pasos, prob;
    int norte, sur, este, oeste, contador;
    contador = 0;
    //pido las dimensiones de la matriz
    cout << "Introduzca las dimensiones de la matriz s: ";
    cin >> N;

    cout << endl;

    //pido la temperatura
    cout << "Introduzca un valor de T (entre 0.0 y 5.0): ";
    cin >> T;

    cout << "Introduzca los pasos Montecarlo: ";
    cin >> pasos;

    //calculo beta
    beta = 1.0/T;

    //relleno la matriz (estado inicial) todos los espines igual a 1
    for(i=1;i<=N;i++)
    {
        for(j=1;j<=N;j++)
        {
            s[i][j]=1;
        }
    }

    double M = N*N*pasos;
    int k = 0;

    //hago N^2 intentos de cambio de Spin (la cantidad de spines en la matriz)
    while(k<M)
    {
        //escojo una posicion  aleatoria de la matriz
        n = rand() % N + 1; //en x
        m = rand() % N + 1; //en y
   // cout << n << " " << m << endl;

        //chequeamos condiciones de contorno periodicas para n y m
        norte = m + 1;
        if(norte > N)
        {
            norte = 1;
        }

        sur = m - 1;
        if(sur < 1)
        {
            sur = N;
        }

        este = n + 1;
        if(este > N)
        {
            este = 1;
        }

        oeste = n - 1;
        if(oeste < 1)
        {
            oeste = N;
        }

        //obtengo dE que es la energía del estado elegido al azar
        dE = 2*s[n][m]*(s[este][m]+s[oeste][m]+s[n][norte]+s[n][sur]);
        //cout << "dE " << dE << endl;

        //obtengo la probabilidad de ese estado
        prob = exp(-dE*beta);
        P=minimo(1.0,prob);
        //cout << "P " << P << endl;

        //genero numero aleatorio hi, entre cero y 1 y comparo con P
        hi = (double)rand()/RAND_MAX;

        //cout << "chi " <<  hi << endl;

        if(hi < P)
        {
            s[n][m] = -s[n][m];
           // MuestraMatriz(s,N);
            contador++;
        }

        k++;
    }

    cout << contador << " " << k << endl;
    //calculo la magnetizacion media
    cout << magnetizacion(N,s) << endl;

    //despliego la matriz por pantalla
    MuestraMatriz(s,N);

    cout << endl;

   // cout << "Energia = " << E << endl;

    return 0;
}

void MuestraMatriz(int s[100][100], int N)
{
    int i,j;

    for(i=1;i<=N;i++)
    {
        for(j=1;j<=N;j++)
        {
            cout << s[i][j] << " ";
        }
        cout << endl;
    }

    return;
}

double Energia(int s[100][100], int N)
{
    int i,j;
    double E = 0.0;

    for(i=1;i<=N;i++)
    {
        for(j=1;j<=N;j++)
        {
            E+=s[i][j]*(s[i][j+1]+s[i][j-1]+s[i+1][j]+s[i-1][j]);
        }
    }

    return -E*0.5;
}

double minimo(double a, double b)
{
    double minimo;

    if(a>b)
    {
        minimo = b;
    }
    else
    {
        minimo = a;
    }

    return minimo;
}

double magnetizacion(int N, int s[100][100])
{
    double magnetizacion = 0.0;
    int i,j;

    for(i=1;i<=N;i++)
    {
        for(j=1;j<=N;j++)
        {
            magnetizacion += s[i][j];
        }
    }

    return abs(magnetizacion/(N*N));
}
