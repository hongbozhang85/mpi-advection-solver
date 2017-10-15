//parallel 2D advection solver module
// written for COMP4300/8300 Assignment 1 
// v1.0 15 Mar  

extern int M_loc, N_loc; // local advection field size (excluding halo) 
extern int M0, N0;       // local field element (0,0) is global element (M0,N0)
extern int P0, Q0;      // 2D process id (P0, Q0) in P x Q process grid 

//sets up parallel parameters above
void initParParams(int M, int N, int P, int Q, int verbosity);

// evolve advection over r timesteps, with (u,ldu) storing the local filed
void parAdvect(int r, double *u, int ldu);

// overlap communication variant
void parAdvectOverlap(int r, double *u, int ldu);

// wide halo variant
void parAdvectWide(int r, int w, double *u, int ldu);

// extra optimization variant
void parAdvectExtra(int r, double *u, int ldu);
