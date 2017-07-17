#include <cmath>

using namespace std;
// funcoes para um lattice D2Q9
/*
0\3/6
1-4-7
2/5\8
*/

double cs2 = 1.0 /3 ;
int cx[9] = {-1,-1,-1,   0 , 0, 0   ,1 ,1, 1} ;
int cy[9] = { 1, 0,-1,   1 , 0,-1   ,1, 0,-1} ;
double w[9] = { 1./36, 1./9,  1./36,     1./9,   4./9,    1./9,       1./36,  1./9,   1./36};


////////////////////////////////////////// Operacoes locais //////////////////////////////////

double calcDensidade(vector<double> &f){
double rho = 0;
for(int i =0  ; i<9 ; i++ ){rho = rho + f[i];}
return rho;
};

double calcVelX(vector<double> &f, double rho){
double Ux;
Ux= ( f[6]+f[7]+f[8] - f[0] -f[1] - f[2]  )/rho;
return Ux;
};

double calcVelY(vector<double> &f, double rho){
double Uy;
Uy = (  f[0]+f[3]+f[6] -f[2] -f[5] -f[8]   )/rho;
return Uy;
};

double calcEquilibrio(double rho, double Ux, double Uy, int q){
double a = (cx[q]*Ux+ cy[q]*Uy ) / cs2; //
double b = (Ux*Ux+Uy*Uy) / cs2;
double f_eq;

f_eq = w[q] * rho* ( 1 + a + a*a /2 - b/2);

return f_eq;
};
////////////////////// Esquemas de Pressao ////////////////////////////////////
/*
void PressaoZuo(vector<double> &f, double rho){
double ux;
ux = 1 - ( (f[4]+f[3]+f[5])+2*(f[1]+f[0]+f[2]) ) / rho;
f[6] = f[2] - 1/2 * (f[3]-f[5]) + 1/6 * rho * ux;
f[7] = f[1] + 2/3 * rho * ux;
f[8] = f[0] + 1/2 * (f[3]-f[5]) + 1/6 * rho * ux;
}
*/
double PressaoLuo(double rhoIn, double rho0, double Ux, double Uy, int q){
double cl = 1.0  ;

double feq;
double a = (cx[q]*Ux+ cy[q]*Uy ) / cl;
double b = (Ux*Ux+Uy*Uy) / (cl*cl);

feq = w[q]*(rhoIn + rho0*(3*a+9/2*a*a-3/2*b));
return feq;
}

/////////////////////////////////// Operacoes no lattice //////////////////////////////////////

void DefineLattice(vector<vector<vector<double> > > &lattice,int nx,int ny){
lattice.resize(ny);
  for (int j = 0; j < ny; j++){
    lattice[j].resize(nx);
        for (int i= 0; i < nx; ++i){
            lattice[j][i].resize(9);}};
};

void IniciaLatticeEq(vector<vector<vector<double> > > &lattice,int x0, int x1, int y0, int y1, double rho, double Ux, double Uy){
for (int j = y0; j< y1; j++){
        for (int i = x0; i<x1; i++){
            for(int q = 0; q < 9; q++){
         lattice[j][i][q] = calcEquilibrio(rho, Ux, Uy,q);
            }}}
};

void DefineBarreiras(vector<vector<bool> >  &Barreiras,int nx,int ny){
Barreiras.resize(ny);
for(int j = 0; j<ny; j++){
   Barreiras[j].resize(nx);
    };
for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        Barreiras[j][i] = false;
    }}
};






