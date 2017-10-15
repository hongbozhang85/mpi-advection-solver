// serial 2D advection solver module
// written by Peter Strazdins, Mar 17 for COMP4300/8300 Assignment 1 
// v1.0 15 Mar 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> // sin(), fabs()

#include "serAdvect.h"


// advection parameters
static const double CFL = 0.25;   // CFL condition number
static const double Velx = 1.0, Vely = 1.0; //advection velocity
static double dt;                 //time for 1 step
static double deltax, deltay;     //grid spacing//

void initAdvectParams(int M, int N) {
  assert (M > 0 && N > 0); // advection not defined for empty grids
  deltax = 1.0 / N;
  deltay = 1.0 / M;
  dt = CFL * (deltax < deltay? deltax: deltay);
  // printf("dx=%e dy=%e dt=%e\n", deltax, deltay, dt);
}

static double initCond(double x, double y, double t) {
  x = x - Velx*t;
  y = y - Vely*t;
  return (sin(4.0*M_PI*x) * sin(2.0*M_PI*y)) ;
}

void initAdvectField(int M0, int N0, int m, int n, double *u, int ldu) {
  int i, j;
  for (j=0; j < n; j++) {
    double x = deltax * (j + N0);
    for (i=0; i < m; i++) {
      double y = deltay * (i + M0);
      V(u, i, j) = initCond(x, y, 0.0);
    }
  }
}

double errAdvectField(int r, int M0, int N0, int m, int n, double *u, int ldu){
  int i, j;
  double err = 0.0;
  double t = r * dt;
  for (j=0; j < n; j++) {
    double x = deltax * (j + N0);
    for (i=0; i < m; i++) {
      double y = deltay * (i + M0);
      err += fabs(V(u, i, j) - initCond(x, y, t));
    }
  }
  return (err);
}

const int AdvFLOPsPerElt = 23; //count 'em

// uses the Lax-Wendroff method
void updateAdvectField(int m, int n, double *u, int ldu) {
  int i, j;
  double sx = 0.5 * Velx / deltax, sy = 0.5 * Vely / deltay;
  double *ut; int ldut = m+1;
  ut = malloc(ldut*(n+1)*sizeof(double));
  assert(ut != NULL);
  for (j=0; j < n+1; j++)     
    for (i=0; i < m+1; i++) {
      V(ut,i,j) = 0.25*(V(u,i,j) + V(u,i-1,j) + V(u,i,j-1) + V(u,i-1,j-1))
	   -0.5*dt*(sy*(V(u,i,j) + V(u,i,j-1) - V(u,i-1,j) - V(u,i-1,j-1)) +
	 	    sx*(V(u,i,j) + V(u,i-1,j) - V(u,i,j-1) - V(u,i-1,j-1)));
    }
  double dtdx = 0.5 * dt / deltax, dtdy = 0.5 * dt / deltay;
  for (j=0; j < n; j++)     
    for (i=0; i < m; i++) {
      V(u, i, j) +=  
	- dtdy * (V(ut,i+1,j+1) + V(ut,i+1,j) - V(ut,i,j) - V(ut,i,j+1))
	- dtdx * (V(ut,i+1,j+1) + V(ut,i,j+1) - V(ut,i,j) - V(ut,i+1,j));
    }
  free(ut);
} //updateAdvectField()


void printAdvectField(int rank, char *label, int m, int n, double *u, int ldu){
  int i, j;
  printf("%d: %s\n", rank, label);
  for (i=0; i < m; i++) {
    printf("%d: ", rank);  
    for (j=0; j < n; j++) 
      printf(" %+0.2f", V(u, i, j));
    printf("\n");
  }
}
