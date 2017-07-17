/*
Fluid flow and acoustics using the Lattice Boltzmann Method
air flow through glottis
Author:Wagner Rupp
*/
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

#include "funcoesD2Q9.h"
#include "Geometrias.h"
#include "Conversor.h"

using std::vector;


int main(){

const int nx = 467;
const int ny = 561;

int tmax = 100000;
const double Omega = 1.9898; // (1/tau)

// Dados do fluido ///

const double C_r = 330;         // [m/s]
const double rho_r0 = 1.2;            // [kg/m3]
const double visc_Kin = 1.5e-5;         // [m2/s]

const double DeltaPressao = 294.5;  // [Pa]

// Parametros D2Q9

vector<double> f; f.resize(9);
const double C_l = 1./sqrt(3);

int cx[9] = {-1,-1,-1,   0 , 0, 0   ,1 ,1, 1} ;
int cy[9] = { 1, 0,-1,   1 , 0,-1   ,1, 0,-1} ;
double w[9] = { 1./36, 1./9,  1./36,     1./9,   4./9,    1./9,       1./36,  1./9,   1./36};
double f_NE, f_0;
double rho, Ux, Uy, U, f_eq;

// Densidade e pressao

double rho_l0 = 1;                    // [lattice]
double dp = rho_r0 / rho_l0 ;         // [kg/m3]
double P_r0 = C_r* C_r * rho_r0;      // [Pa]
double rho_r1 = ( P_r0 + DeltaPressao )/( C_r * C_r );

double rhoIn =  rho_r1/dp;            // [lattice]
double rhoOut =  1 ;                  // [lattice]

const double dx = 3.077e-5;           // [m]
const double dt = C_l/C_r * dx ;      // [s]

// Condicoes iniciais
double rho0 = rho_l0;
double Ux0 = 0;
double Uy0 = 0;

// Iteradores
int j , i, q;
int indexJ;
int indexI;

// Saida de dados
int display = 1;
ofstream arquivoD[10];
ofstream arquivoU[10];
ofstream arquivoB;
std::string nomeDensi ="Densidade_";
std::string nomeVelo = "Velocidade_";
std::string strD,strU;

//Alocando o lattice
std::cout << "Alocando memoria "<<std::endl;

vector<vector<vector<double> > > lattice;
vector<vector<vector<double> > > f_strm;
vector<vector<bool> > Barreiras;

DefineLattice(lattice,nx,ny);
DefineLattice(f_strm,nx,ny);
DefineBarreiras( Barreiras, nx, ny);

///////////////////////  Condicoes Iniciais  //////////////////////////
std::cout << "Implementando Condicoes iniciais "<<std::endl;
IniciaLatticeEq(lattice,0, nx, 0, ny, rho0, Ux0, Uy0);
//IniciaLatticeEq(lattice,0,  2, 1, ny, rhoIn, Ux0, Uy0);
////////////////////////////// Barreiras //////////////////////////////

std::cout << "Implementando Barreiras "<<std::endl;

PregaConv(dx, nx, ny, Barreiras );  // Conv, Uni ou Div

arquivoB.open("Barreiras.csv");
for( j = 0; j< ny; j++){
        for( i=0; i<nx ; i++){
                arquivoB << Barreiras[j][i] << " ; ";
        }
        arquivoB<<std::endl;
}
arquivoB.close();


//////////////////////// Loop Principal //////////////////////////////
std::cout << "Iniciando Loop principal "<<std::endl;

for(int t = 0;  t<=tmax ; t++){
                    //////////// Colisao ///////////
for( j = 0; j<ny; j++){
for( i=0; i<nx ; i++){
    if(Barreiras[j][i]==false ){
        for( q=0; q<9; q++){f[q] = lattice[j][i][q];}
        rho = calcDensidade(f);
        Ux = calcVelX(f,rho);
        Uy = calcVelY(f,rho);
        for( q=0; q<9 ; q++){
            f_eq = calcEquilibrio(rho,Ux,Uy,q);
            f_strm[j][i][q] = - Omega*(f[q] - f_eq) + f[q];}

    ///////////Pressao//////////////

            if(i==0 && j>1 && j<ny-2 ){
              for( q=0; q<9; q++){f[q] = lattice[j][i+1][q];}
                rho = calcDensidade(f);
                Ux = calcVelX(f,rho);
                Uy = calcVelY(f,rho);
              for(q=0; q<9; q++){
                f_0 = PressaoLuo(rhoIn,rho0, Ux , Uy, q);
                f_NE = f[q] - calcEquilibrio(rho,Ux,Uy,q);
                f_strm[j][i][q] = f_0 + (1 - Omega)*f_NE;
                        }}


            if(i==nx-1 && j>1 && j<ny-2){
              for( q=0; q<9; q++){f[q] = lattice[j][i-1][q];}
                rho = calcDensidade(f);
                Ux = calcVelX(f,rho);
                Uy = calcVelY(f,rho);
               for(q=0; q<9; q++){
                f_NE = f[q] - calcEquilibrio(rho,Ux,Uy,q);
                f_0 = PressaoLuo(rhoOut, rho0, Ux , Uy, q);
                f_strm[j][i][q] = f_0 + (1 - Omega)*f_NE;
                        }}


                      }
    else{
      for( q=0; q<9; q++){f_strm[j][i][q] = lattice[j][i][q];}
        }


 }
 }

                ////////// Propagacao ///////////

    for(j=1;j<ny-1;j++){
        for(i=1; i<nx-1; i++){
            for(q=0;q<9; q++){

              indexJ = j - cy[q];
              indexI = i - cx[q];

 lattice[j][i][q]  = f_strm[indexJ][indexI][q];




                              }


                //////// Retorno Barreiras //////////////


  if(i==0){
    lattice[j][i][0] = f_strm[j][i][8];
    lattice[j][i][1] = f_strm[j][i][7];
    lattice[j][i][2] = f_strm[j][i][6];
                    }

  if(j==1){
    lattice[j][i][8] = f_strm[j][i][0];
    lattice[j][i][5] = f_strm[j][i][3];
    lattice[j][i][2] = f_strm[j][i][6];
                    }

  if( j==ny-2 ){
    lattice[j][i][0] = f_strm[j][i][8];
    lattice[j][i][3] = f_strm[j][i][5];
    lattice[j][i][6] = f_strm[j][i][2];
                    }

   if(Barreiras[j][i]==true ){
    lattice[j][i][0] = f_strm[j][i][8];
    lattice[j][i][1] = f_strm[j][i][7];
    lattice[j][i][2] = f_strm[j][i][6];
    lattice[j][i][3] = f_strm[j][i][5];

    lattice[j][i][5] = f_strm[j][i][3];
    lattice[j][i][6] = f_strm[j][i][2];
    lattice[j][i][7] = f_strm[j][i][1];
    lattice[j][i][8] = f_strm[j][i][0];
                    }




        }
    }

            ///////////////// Output /////////////////

    if(t==tmax/10*display){
            std::cout << display<< 0 <<"% concluido"<< std::endl;
            strD=nomeDensi+ std::to_string(t)+".csv"; // Mudar esse char(t)
            arquivoD[display-1].open(strD.c_str());
            strU=nomeVelo+ std::to_string(t)+".csv";
            arquivoU[display-1].open(strU.c_str());

                for (j=0;j<ny;j++){
                        for(i=0;i<nx;i++){
                            for(q=0;q<9;q++){f[q]=lattice[j][i][q];}

                            rho = calcDensidade(f);
                            Ux = calcVelX(f,rho);
                            Uy = calcVelY(f,rho);
                            U = sqrt(Ux*Ux+Uy*Uy);

                            arquivoD[display-1] << rho << " ; " ;
                            arquivoU[display-1] << U << " ; " ;
                                        }
                        arquivoD[display-1] << endl;
                        arquivoU[display-1] << endl;
                                }
                arquivoD[display-1].close();
                arquivoU[display-1].close();
                display++;
                }

}

return 0;
};
