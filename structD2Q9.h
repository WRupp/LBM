#include <cmath>

struct d2q9 {
int nx;
int ny;

double c_l = 1./sqrt(3);

double Omega;

int cx[9] = {-1,-1,-1,   0 , 0, 0   ,1 ,1, 1} ;
int cy[9] = { 1, 0,-1,   1 , 0,-1   ,1, 0,-1} ;
double w[9] = { 1./36, 1./9,  1./36,     1./9,   4./9,    1./9,       1./36,  1./9,   1./36};
};
