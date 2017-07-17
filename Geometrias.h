#include <vector>
#include <cmath>

using namespace std;

// funcoes para geracao de barreiras
/*
0\3/6
1-4-7
2/5\8
*/

void criaCirculo(int Cx,int Cy, double R, int nx, int ny, vector<vector<bool> > &Barreiras){
int circulo;
double R2;
R2 = R*R;

for(int j=0;j<ny;j++){
	for(int i=0;i<nx; i++){
	    circulo = (i-Cx)*(i-Cx)+(j-Cy)*(j-Cy);
		if( circulo<R2){Barreiras[j][i] = true;}}}
};

void criaLinha(double x1, double y1, double x2, double y2, vector<vector<bool> > &Barreiras){
// Uma funcao simples de linha, existe o problema das "frestas", uma opcao seria implementar
// o algoritmo de linhas de Bresenham

double a,b;
int index1,index2,y3;
if(x1!=x2){
    a = (y2-y1)/(x2-x1);
    b = y2-a*x2;
        if(y2>=y1){index1=y1; index2=y2;}
        else{index1=y2; index2=y1;}

     for(int j=index1; j<=index2; j++){
        for(int i=x1;i<=x2; i++){
            y3 = a*i+b;
            if(j == floor(y3) || j==ceil(y3)) {
                Barreiras[j][i]=true;
                }}}}
else{  for(int j=y1;j<=y2;j++){Barreiras[j][x1]=true;}  }
};

void linhaB(int x0, int y0, int x1, int y1, vector<vector<bool> > &Barreiras){
// o Algoritmo de Linhas de Bresenham, adaptado

double deltaX = x1-x0;
double deltaY = y1-y0;
double ErroDelta = abs(deltaY/deltaX );
double erro = ErroDelta-0.5;
int y = y0;
int ymax = Barreiras.size()-2;
int mid = Barreiras.size()/2-1;

if(x1!=x0){
    for(int x=x0;x<=x1;x++){
        y = (deltaY/deltaX)*(x-x0)+y0;
        Barreiras[y][x]=true;
        Barreiras[y][x+1]=true;
        erro = erro + ErroDelta;
            if(erro>=1 && y <ymax && y>=mid ){
                y=y+1;
                Barreiras[y][x]=true;
                Barreiras[y][x+1]=true;
                erro = erro - 1.0;
                        }
            if(erro>=1 && y > 1 && y<mid  ){
                y=y-1;
                Barreiras[y][x]=true;
                Barreiras[y][x+1]=true;
                erro = erro - 1.0;
                        }

                            }
          }
else{  for(int j=y0;j<=y1;j++){Barreiras[j][x1]=true;}  }
};

////////////////////////// PREGAS /////////////////////////////////////////

void PregaConv(double dx,int nx, int ny, vector<vector<bool> > &Barreiras ){
// dx = xr /xl ; xl=xr/dx;

int Ox,Oy;
int Cx1,Cy1,R1;
int Cx2,Cy2,R2;
int Ax1,Ay1,Ax2,Ay2;
int Bx1,By1,Bx2,By2;
int Dx1,Dy1,Dx2,Dy2;
int Ex,Ey1,Ey2;

int cy1,cy2;
int ay1,ay2;
int by1,by2;
int dy1,dy2;
int ey1,ey2;

Ox=ceil(9.6e-3/dx)-1;
Oy=ceil(ny/2-1)-1;

R1 = ceil(1.5e-3/dx);
Cx1 = Ox + 0;
Cy1 = Oy - ceil(1.89e-3/dx) ;
cy1 = Oy + ceil(1.89e-3/dx) ;

R2 = ceil(0.908e-3/dx);
Cx2 = Ox + ceil(2.093e-3/dx);
Cy2 = Oy - ceil(1.108e-3/dx);
cy2 = Oy + ceil(1.108e-3/dx);

Ax1 = ceil(2.0e-3/dx); // REF GLOBAL - TROCAR PRA LOCAL
Ax2 = Ox - ceil(1.059e-3/dx);

Ay1 = 0; // REF GLOBAL - TROCAR PRA LOCAL
ay1 = ny -1;
Ay2 = Oy - ceil(0.823e-3/dx);
ay2 = Oy + ceil(0.823e-3/dx);

Bx1 = Ox + ceil(0.131e-3/dx);
Bx2 = Ox + ceil(2.013e-3/dx);

By1 = Oy - ceil(0.391e-3/dx)-1;
by1 = Oy + ceil(0.391e-3/dx)+1;

By2 = Oy - ceil(0.203e-3/dx)-1;
by2 = Oy + ceil(0.203e-3/dx)+1;

Ex  = Cx2+ R2;
Ey1 = ceil(0.8e-3/dx);
ey1 = cy2;
Ey2 = Cy2;
ey2 = ny - ceil(0.8e-3/dx)-1 ;

Dx1 = Cx2+R2;
Dx2 = Cx2+R2+ceil(0.8e-3/dx)+1;
Dy1 = ceil(0.8e-3/dx); // REF GLOBAL - TROCAR PRA LOCAL
dy1 = ny - ceil(0.8e-3/dx)-1; // REF GLOBAL - TROCAR PRA LOCAL
Dy2 = 0; // REF GLOBAL - TROCAR PRA LOCAL
dy2 = ny - 1; // REF GLOBAL - TROCAR PRA LOCAL

criaCirculo(Cx1,Cy1,R1,nx,ny,Barreiras);
criaCirculo(Cx2,Cy2,R2,nx,ny,Barreiras);

linhaB(Ax1,Ay1,Ax2,Ay2, Barreiras);
linhaB(Bx1,By1,Bx2,By2, Barreiras);
linhaB(Dx1,Dy1,Dx2,Dy2, Barreiras);
linhaB(Ex,Ey1,Ex,Ey2, Barreiras);

criaCirculo(Cx1,cy1,R1,nx,ny,Barreiras);
criaCirculo(Cx2,cy2,R2,nx,ny,Barreiras);

linhaB(Ax1,ay1,Ax2,ay2, Barreiras);
linhaB(Bx1,by1,Bx2,by2, Barreiras);
linhaB(Dx1,dy1,Dx2,dy2, Barreiras);
linhaB(Ex,ey1,Ex,ey2, Barreiras);
};

void PregaUni(double dx,int nx, int ny, vector<vector<bool> > &Barreiras ){
// dx = xr /xl ; xl=xr/dx;

int Ox,Oy;
int Cx1,Cy1,R1;
int Cx2,Cy2,R2;
int Ax1,Ay1,Ax2,Ay2;
int Bx1,By1,Bx2,By2;
int Dx1,Dy1,Dx2,Dy2;
int Ex,Ey1,Ey2;

int cy1,cy2;
int ay1,ay2;
int by1,by2;
int dy1,dy2;
int ey1,ey2;

Ox=ceil(9.6e-3/dx)-1;
Oy=ceil(ny/2-1)-1;

R1 = ceil(1.5e-3/dx);
Cx1 = Ox + 0;
Cy1 = Oy - ceil(1.7e-3/dx) ;
cy1 = Oy + ceil(1.7e-3/dx) ;

R2 = ceil(0.987e-3/dx);
Cx2 = Ox + ceil(2.013e-3/dx);
Cy2 = Oy - ceil(1.18e-3/dx);
cy2 = Oy + ceil(1.18e-3/dx) ;

Ax1 = ceil(2.0e-3/dx);    // REF GLOBAL - TROCAR PRA LOCAL
Ax2 = Ox - ceil(1.149e-3/dx);
Ay1 = 0;
Ay2 = Oy - ceil(0.736e-3/dx);
ay1 = ny -1;
ay2 = Oy + ceil(0.736e-3/dx);

Bx1 = Ox + 0;
Bx2 = Ox + ceil(2.013e-3/dx);
By1 = Oy - ceil(0.2e-3/dx);
by1 = Oy + ceil(0.2e-3/dx);
By2 = Oy - ceil(0.2e-3/dx);
by2 = Oy + ceil(0.2e-3/dx);

Ex  = Cx2+R2;
Ey1 = ceil(0.8e-3/dx); // REF GLOBAL - TROCAR PRA LOCAL
ey1 = cy2 -1;
Ey2 = Cy2;
ey2 = ny - ceil(0.8e-3/dx)-1 ; // REF GLOBAL - TROCAR PRA LOCAL



Dx1 = Cx2+R2;
Dx2 = Cx2+R2+ceil(0.8e-3/dx)+1;
Dy1 = ceil(0.8e-3/dx); // REF GLOBAL - TROCAR PRA LOCAL
dy2 = ny - 1; // REF GLOBAL - TROCAR PRA LOCAL
Dy2 = 0;// REF GLOBAL - TROCAR PRA LOCAL
dy1 = ny - ceil(0.8e-3/dx)-1; // REF GLOBAL - TROCAR PRA LOCAL


criaCirculo(Cx1,Cy1,R1,nx,ny,Barreiras);
criaCirculo(Cx2,Cy2,R2,nx,ny,Barreiras);

criaCirculo(Cx1,cy1,R1,nx,ny,Barreiras);
criaCirculo(Cx2,cy2,R2,nx,ny,Barreiras);

linhaB(Ax1,Ay1,Ax2,Ay2, Barreiras);
linhaB(Bx1,By1,Bx2,By2, Barreiras);
linhaB(Dx1,Dy1,Dx2,Dy2, Barreiras);
linhaB(Ex,Ey1,Ex,Ey2, Barreiras);

linhaB(Ax1,ay1,Ax2,ay2, Barreiras);
linhaB(Bx1,by1,Bx2,by2, Barreiras);
linhaB(Dx1,dy1,Dx2,dy2, Barreiras);
linhaB(Ex,ey1,Ex,ey2, Barreiras);
}

void PregaDiv(double dx,int nx, int ny, vector<vector<bool> > &Barreiras ){

int Ox,Oy;
int Cx1,Cy1,R1;
int Cx2,Cy2,R2;
int Ax1,Ay1,Ax2,Ay2;
int Bx1,By1,Bx2,By2;
int Dx1,Dy1,Dx2,Dy2;
int Ex,Ey1,Ey2;

int cy1,cy2;
int ay1,ay2;
int by1,by2;
int dy1,dy2;
int ey1,ey2;

Ox=ceil(9.6e-3/dx)-1;
Oy=ceil(ny/2-1)-1;

R1 = ceil(1.5e-3/dx);
Cx1 = Ox + 0;
Cy1 = Oy - ceil(1.7e-3/dx) ;
cy1 = Oy + ceil(1.7e-3/dx) ;

R2 = ceil(1.08e-3/dx);
Cx2 = Ox + ceil(1.92e-3/dx);
Cy2 = Oy - ceil(1.446e-3/dx);
cy2 = Oy + ceil(1.446e-3/dx) ;

Ax1 = ceil(2.0e-3/dx);    // REF GLOBAL - TROCAR PRA LOCAL
Ax2 = Ox - ceil(1.149e-3/dx);
Ay1 = 0;
Ay2 = Oy - ceil(0.736e-3/dx);
ay1 = ny -1;
ay2 = Oy + ceil(0.736e-3/dx);

Bx1 = Ox + 0;
Bx2 = Ox + ceil(2.014e-3/dx);
By1 = Oy - ceil(0.206e-3/dx);
by1 = Oy + ceil(0.206e-3/dx);
By2 = Oy - ceil(0.370e-3/dx);
by2 = Oy + ceil(0.370e-3/dx);

Ex  = Cx2+R2;
Ey1 = ceil(0.8e-3/dx); // REF GLOBAL - TROCAR PRA LOCAL
ey2 = ny - ceil(0.8e-3/dx)-1 ; // REF GLOBAL - TROCAR PRA LOCAL
Ey2 = Cy2;
ey1 = cy2 -1;

Dx1 = Cx2+R2;
Dx2 = Cx2+R2+ceil(0.8e-3/dx)+1;
Dy1 = ceil(0.8e-3/dx); // REF GLOBAL - TROCAR PRA LOCAL
dy2 = ny - 1; // REF GLOBAL - TROCAR PRA LOCAL
Dy2 = 0;// REF GLOBAL - TROCAR PRA LOCAL
dy1 = ny - ceil(0.8e-3/dx)-1; // REF GLOBAL - TROCAR PRA LOCAL


criaCirculo(Cx1,Cy1,R1,nx,ny,Barreiras);
criaCirculo(Cx2,Cy2,R2,nx,ny,Barreiras);

criaCirculo(Cx1,cy1,R1,nx,ny,Barreiras);
criaCirculo(Cx2,cy2,R2,nx,ny,Barreiras);

linhaB(Ax1,Ay1,Ax2,Ay2, Barreiras);
linhaB(Bx1,By1,Bx2,By2, Barreiras);
linhaB(Dx1,Dy1,Dx2,Dy2, Barreiras);
linhaB(Ex,Ey1,Ex,Ey2, Barreiras);

linhaB(Ax1,ay1,Ax2,ay2, Barreiras);
linhaB(Bx1,by1,Bx2,by2, Barreiras);
linhaB(Dx1,dy1,Dx2,dy2, Barreiras);
linhaB(Ex,ey1,Ex,ey2, Barreiras);
};

