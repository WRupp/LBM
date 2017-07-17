#ifndef CONVERSOR_H_INCLUDED
#define CONVERSOR_H_INCLUDED

#include <cmath>

#endif // CONVERSOR_H_INCLUDED

///////////////////////////////////// Conversão entre sistemas ////////////////////////////////////

int x_r2l(double dx, double distancia){
return ceil(distancia/dx) ;
};

int x_l2r(double dx,int x_lat){
return x_lat*dx;
};

double visc2Tau(double Visc_real,double dx, double C_real){
double tau;
double visc_lat;
double C_l = 1./sqrt(3);

visc_lat = Visc_real * (C_l/C_real) * (1.0/dx);
// Visc_kin = cs2(tau-1/2)
tau = (C_l*C_l)*visc_lat + 0.5;
return tau;
}

//Pressao real->lattice
double tau2visc(double tau,double dx, double C_real){
double visc_lat;
double visc_real;
double C_l = 1./sqrt(3);

visc_lat = (tau-0.5)/3;
visc_real = visc_lat*dx*(C_real/C_l);
return visc_real;
}

//Pressao real->lattice
double P_r2l(double P_r,double rho_r,double C_r){
return P_r/rho_r * 1/C_r;
}


//Converte velocidade de lattice p/ real
//Mach = U /C;
double U_l2r(double Ui , double C_real){
double U_real; double Cl = 1./sqrt(3);
U_real = (C_real / Cl) * Ui;
return U_real;
};

//Converte velocidade real p/lattice
double U_r2l(double U_real , double C_real){
double U_lat; double Cl = 1./sqrt(3);
U_lat = (Cl / C_real) * U_real;
return U_lat;
};

// Re = L*U / Visc_kin
double calcReynolds(double L, double U, double Visc_kin){
double Re;
Re = L*U/Visc_kin;
return Re;
}

