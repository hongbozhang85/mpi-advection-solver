// serial 2D advection solver module
// written by Peter Strazdins, Mar 17 for COMP4300/8300 Assignment 1 
// v1.0 15 Mar

// number of FLOPs to update a single element in the advection function 
extern const int AdvFLOPsPerElt; 

// initializes the advection paamters for a global M x N field 
void initAdvectParams(int M, int N);

// access element (i,j) of array u with leading dimension ldu
#define V(u, i, j) u[(i) + (ld##u)*(j)]

//update 1 timestep for the local advection, without updating halos
//  the m x n column major array (u,ldu) is updated.
//  Assumes a halo of width 1 are around this array;
//  the corners of the halo are at a[-1,-1], a[-1,n], a[m,-1] and a[m,n]
void updateAdvectField(int m, int n, double *u, int ldu);

// initialize (non-halo elements) of a m x n local advection field (u, ldu)
//    local element [0,0] is element [M0,N0] in the global field
void initAdvectField(int M0, int N0, int m, int n, double *u, int ldu);

// sum errors in an m x n local advection field (u, ldu) after r timesteps 
//    local element [0,0] is element [M0,N0] in the global field 
double errAdvectField(int r, int M0, int N0, int m, int n, double *u, int ldu);

// print out  m x n local advection field (u, ldu) 
void printAdvectField(int rank, char *label, int m, int n, double *u, int ldu);

