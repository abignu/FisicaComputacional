#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#define PI 3.1415926535

void CalcularAceleraciones(double r[10][2], double m[10], double a[10][2], int N);
double Hamiltoniano(double m[10], double v[10][2], double r[10][2], int N);
void CalculaPeriodo(double r[10][2], double v[10][2], int N, double v_angular, double periodo[10]);

using namespace std;

int main(void)
{
    //Declaro variables, r[1][1] posicion x r[1][2] posicion y de un planeta. Cada planeta tendra una posicion que se irá renovando con los bucles.

    int i,j,k,N;
    double Ms,G,r[10][2],v[10][2],a[10][2],h,t,c,t_aux,l_angular[10],t_,w[10][2], M, p[10], v_angular, periodo[10], H[10];
    double Tmax, Msistema, cm[2];
    ofstream planetas("planetas.txt");
    ifstream posiciones("posiciones.txt");
    ifstream velocidades("velocidades.txt");
    ofstream mercurio("mercurio.txt"), venus("venus.txt"), tierra("tierra.txt"), jupiter("jupiter.txt"), saturno("saturno.txt"), urano("urano.txt"), neptuno("neptuno.txt"), pluton("pluton.txt"), marte("marte.txt");
    ofstream sol("sol.txt"), mlineal("m_lineal.txt"), energia("energia.txt");

    //masa del sol y constantes universal
    Ms=1.99e30;
    G=6.67e-11;
    M=1.0;
    Msistema=0.0;

    //constante de distancia en unidades astronomicas. Todas las distancias estarán en unidades c. 4c serian 4 unidades astronomicas
    c=1.496e11;

    //t_ es la variable para transformar las coordenadas de tiempo
    t_aux=sqrt(((G)*(Ms))/(pow(c,3.0)));

    //h es el paso medido en unidades de tiempo, ex[i] es la excentricidad del planeta i, igual que m[i] la masa.
    //el tiempo ira en años terrestres, la distancia al sol en unidades astronómicas y la masa en masas solares.
    //introducimos la masa de cada planeta, ya normalizada, y la excentricidad.

    double m[] = {1.0, 1.666e-7,	2.459e-6, 3.015e-6,	3.240e-7, 9.590e-4, 2.86e-4, 4.384e-5, 5.151e-5, 6.313e-9};
    double ex[] = {0.0, 0.205, 0.007, 0.017, 0.094, 0.049, 0.057, 0.046,	0.011, 0.244};
    string planeta[] = {"Sol", "Mercurio", "Venus", "Tierra", "Marte", "Júpiter", "Saturno", "Urano", "Neptuno", "Plutón"};

    //pido el numero de cuerpos que interaccionan
    cout << "Introduzca el nº de cuerpos: " << endl;
    cin >> N;

    //definimos las condiciones iniciales. Son 18, dos por planeta. En cuanto al radio, todos arrancarán en el perihelio de su órbita.

    //posiciones componente x,y
    for(i=0;i<=N;i++)
    {
        posiciones >> r[i][0] >> r[i][1];

        //reescalo
        r[i][0]=r[i][0]/c;
        r[i][1]=r[i][1]/c;

    }

    //velocidades componente x,y
    for(j=0;j<=N;j++)
    {
        velocidades >> v[j][0] >> v[j][1];

        //reescalo
        v[j][0]=v[j][0]*sqrt(c/(G*Ms));
        v[j][1]=v[j][1]*sqrt(c/(G*Ms));

    }

    //aquí hacemos el algoritmo de verlet usando r[][], v[][], a[][]
    cout << "Introduzca el tiempo máximo a transcurrir (en segundos): " << endl;
    cin >> Tmax;

    cout << "Introduzca el paso entre planetas (en segundos): " << endl;
    cin >> h;
    t=0.0;
    t_=t*t_aux;

    h=h*t_aux;

    Tmax=Tmax*t_aux;
    //calculo aceleraciones inciales
    CalcularAceleraciones(r,m,a,N);

    // ahora escribo las posiciones a las que les reste el cm
    for(i=0;i<=N;i++)
    {
        if(i==0)
          {sol << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==1)
          {mercurio << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==2)
          {venus << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==3)
          {tierra << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==4)
          {marte << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==5)
          {jupiter << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==6)
          {saturno << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==7)
          {urano << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==8)
          {neptuno << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
        if(i==9)
          {pluton << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << 0 << endl;}
    }
    //ahora el momento lineal inicial
    double momento = 0.0;
    for(i=0;i<=N;i++)
    {
        p[i]=m[i]*sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]);
        momento = momento + p[i];
    }

    //lo escribo en el fichero de momentos
    mlineal << momento << " " << 0 << endl;

    // ahora hago el while del tiempo
    while(t_<Tmax)
    {
        //incremento el tiempo
         t_=t_+h;

        //calculo nuevas r y matriz w con las aceleraciones y las antiguas posiciones y velocidades
        for(i=0;i<=N;i++)
        {
            r[i][0]=r[i][0]+v[i][0]*h+0.5*a[i][0]*h*h;
            r[i][1]=r[i][1]+v[i][1]*h+0.5*a[i][1]*h*h;

            w[i][0]=v[i][0]+h*0.5*a[i][0];
            w[i][1]=v[i][1]+h*0.5*a[i][1];
        }

        //calculo nuevas aceleraciones con las nuevas posiciones
            CalcularAceleraciones(r,m,a,N);

        //calculo nuevas velocidades sabiendo posiciones y nuevas y antiguas aceleraciones

           for(i=0;i<=N;i++)
           {
               v[i][0]=w[i][0]+h*0.5*a[i][0];
               v[i][1]=w[i][1]+h*0.5*a[i][1];

           }

           //las escribo ya restadas junto al periodo

           CalculaPeriodo(r,v,N,v_angular,periodo);

           for(i=0;i<=N;i++)
           {
               if(i==0)
                  {sol << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==1)
                  {mercurio << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==2)
                  {venus << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==3)
                  {tierra << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==4)
                  {marte << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==5)
                  {jupiter << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==6)
                  {saturno << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==7)
                  {urano << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==8)
                  {neptuno << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
                if(i==9)
                  {pluton << r[i][0]-r[0][0] << " " << r[i][1]-r[0][1] << " " << t_ << endl;}
           }

        //verifico el hamiltoniano y las escribo por pantalla
        energia << Hamiltoniano(m,v,r,N) << endl;

        //aca calculo el momento lineal del sistema
        momento = 0.0;
        for(i=0;i<=N;i++)
        {
            //solo componente z
            momento = momento + m[i]*(r[i][0]*v[i][1]-r[i][1]*v[i][0]);
        }


        mlineal << momento << endl;
    }

    return 0;

}

void CalcularAceleraciones(double r[10][2], double m[10], double a[10][2], int N)
{

    int cont1, cont3;
    double sumaX=0.0;
    double sumaY=0.0;

    for(cont1=0;cont1<=N;cont1++)
    {
        sumaX = sumaY = 0.0;
        for(cont3=0;cont3<=N;cont3++)
        {
                if(cont1 != cont3)
                {

                    sumaX=sumaX+(((r[cont1][0]-r[cont3][0])*m[cont3])/pow(pow((r[cont1][0]-r[cont3][0]),2.0)+pow((r[cont1][1]-r[cont3][1]),2.0),1.5));
                    sumaY=sumaY+(((r[cont1][1]-r[cont3][1])*m[cont3])/pow(pow((r[cont1][0]-r[cont3][0]),2.0)+pow((r[cont1][1]-r[cont3][1]),2.0),1.5));

                }
        }
         a[cont1][0]=-sumaX;
         a[cont1][1]=-sumaY;


    }
    return;
}

double Hamiltoniano(double m[10], double v[10][2], double r[10][2], int N)
{
    int k,i;
    double H, V;
    double denominador;
    H=0.0;
    for(k=0;k<=N;k++)
    {
        V=0.0;
        for(i=0;i<=N;i++)
        {
            if((k!=i)&&(i>k))
            {
                denominador=sqrt(pow((r[k][0]-r[i][0]),2.0)+pow((r[k][1]-r[i][1]),2.0));
                V+=(m[k]*m[i])/denominador;
            }
            H+=0.5*m[k]*(v[k][0]*v[k][0]+v[k][1]*v[k][1])-V;
        }
    }

    return H;
}

void CalculaPeriodo(double r[10][2], double v[10][2], int N, double v_angular, double periodo[10])
{
    int k;

    for(k=0;k<=N;k++)
    {
        v_angular=sqrt(v[k][0]*v[k][0]+v[k][1]*v[k][1])/sqrt(r[k][0]*r[k][0]+r[k][1]*r[k][1]);

        periodo[k]=(2*PI)/v_angular;
    }

    return;
}
