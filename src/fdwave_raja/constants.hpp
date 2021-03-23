#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#define POW2(x) ((x)*(x))
#define IDX3(i,j,k) (nz*ny*(i) + nz*(j) + (k))
#define IDX3_l(i,j,k) ((nz+2*lz)*(ny+2*ly)*((i)+lx) + (nz+2*lz)*((j)+ly) + ((k)+lz))
#define IDX3_eta1(i,j,k) ((nz+2)*(ny+2)*((i)+1) + (nz+2)*((j)+1) + ((k)+1))
#define IDX3_eta0(i,j,k) ((nz+2)*(ny+2)*(i) + (nz+2)*(j) + (k))

extern const float _fmax;
extern const float vmin;
extern const float vmax;
extern const float cfl;

typedef long long int llint;
typedef unsigned int uint;

#endif
