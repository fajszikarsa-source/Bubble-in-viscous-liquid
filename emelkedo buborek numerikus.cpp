#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

const double g=-9.8, p0=1e5, M=0.028, R=8.31;
///parameterek
double  V0=1e-8,
        z0=-.6,
        T0=300,
        ro_f=1000,
        eta=.03;

double p(double h){return p0+ro_f*g*h;};
double V(double h){return p(z0)/p(h)*V0;};      ///Boyle-Mariotte tv.
double r(double h){return cbrt(3*V(h)/(4*M_PI));};
double ro(double h){return p(h)*M/(R*T0);};

double egyensulyi_kozelites(double h=z0) {return max(0.0, 2.7*eta/pow(ro_f*g, 2)*pow(4.0*M_PI/(3*p(h)*V(h)), 2.0/3)*(pow(p(h), 5.0/3)-pow(p(-r(h)), 5.0/3)));}

double A(double h, double v){return g*(1.0-ro_f/ro(h))-6.0*M_PI*eta*v*r(h)/(ro(z0)*V0);}    ///gyorsulas hely es sebesseg fuggvenyeben
inline void leptet(double &h, double &v, double &t, const double &dt)      ///masodrendu Runge-Kutta egy lepese
{
    double k11 = dt*v;
    double k21 = dt*A(h, v);
    double k12 = dt*(v+0.5*k21);
    double k22 = dt*A(h+0.5*k11, v+0.5*k21);
    double k13 = dt*(v+0.5*k22);
    double k23 = dt*A(h+0.5*k12, v+0.5*k22);
    double k14 = dt*(v+k23);
    double k24 = dt*A(h+k13, v+k23);
    h += (k11+2*k12+2*k13+k14)/6;
    v += (k21+2*k22+2*k23+k24)/6;
    t+=dt;
}

double sebessegteszt(double dt=egyensulyi_kozelites()*1e-5)     ///szamitogep sebesseget teszteli, hogy tudjunk becslest adni a program futasidejere
{
    double t=0, v=0, h=z0;
    auto start=high_resolution_clock::now();
    for(int i=1e5; --i;) leptet(h, v, t, dt);               ///elkepzelheto, hogy a compiler ezt kioptimalizalja, mert latja, sehol se hasznaljuk a vegeredmenyt
    duration<double> idotartam=high_resolution_clock::now()-start;///ebben az esetben nem fog a futasido elorejelzes jol mukodni
    return idotartam.count();
}

double futasido, Re;
double szimulacio(double dt=egyensulyi_kozelites()*1e-5)
{
    futasido*=1.3;
    bool voltkiiras=false;
    double t=0, v=0, h=z0;
    while(h+r(h)<0){
        leptet(h, v, t, dt);
        if(A(h, v)<0 || v<0) return szimulacio(dt/1.3);         ///ha tul kicsi a lepeskoz ilyen hulyesegek jonnek ki, ilyenkor ujracsinaljuk az egeszet kisebb lepeskozzel
        if(!voltkiiras && t>dt*1000 && futasido>4){             ///valamiert ez a jelenseg neha erthetetlenul erzekeny (kis lepeskozt igenyel): "stiff equation"
            cout<<"\nSzimulacio varhato hatralevo futasideje: "<<(int)futasido<<" masodperc";
            voltkiiras=true;
        }
    }
    Re=ro_f*r(h)*v/eta;
    return t;
}

int main()
{
    futasido=sebessegteszt();
    double t1=egyensulyi_kozelites();
    cout<<"Kozelito keplet eredmenye:\t"<<t1;
    double t2=szimulacio();
    cout<<"\nSzimulacio altal josolt ertek:\t"<<t2<<endl<<setprecision(2)<<abs(t1-t2)/t2*100<<"%-os hiba\n";
    cout<<"Maximalis Reynolds szam a mozgas alatt: "<<fixed<<Re<<'\n';
}
