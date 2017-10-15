// parallel 2D advection solver module
// written for COMP4300/8300 Assignment 1 
// v1.0 15 Mar 

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#include "serAdvect.h"

#define HALO_TAG 100

int M_loc, N_loc; // local advection field size (excluding halo) 
int M0, N0;       // local field element (0,0) is global element (M0,N0)
static int P0, Q0; // 2D process id (P0, Q0) in P x Q process grid 

static int M, N, P, Q; // local store of problem parameters
static int verbosity;
static int rank, nprocs;       // MPI values
static MPI_Comm comm;


//sets up parallel parameters above
void initParParams(int M_, int N_, int P_, int Q_, int verb) {
  M = M_, N = N_; P = P_, Q = Q_;
  verbosity = verb;
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  Q0 = rank % Q;
  N0 = (N / Q) * Q0;
  N_loc = (Q0 < Q-1)? (N / Q): (N - N0); 

  P0 = rank / Q; 
  M0 = (M / P) * P0;
  M_loc = (P0 < P-1)? (M / P): (M - M0); 
} //initParParams()


static void updateBoundary(int w, double *u, int ldu) {
  int i, j;
#  ifndef BLOCK
  MPI_Request reqTR, reqTL, reqFR, reqFL;
#  endif
  // left and right sides of halo
  if (Q == 1) { 
    for (i = w; i < M_loc+w; i++) {
      for (j = 0; j < w; j++) {
        V(u, i, j) = V(u, i, N_loc+j);
        V(u, i, N_loc+w+j) = V(u, i, w+j);
      }
    }
  } else {
    int leftProc = (Q0 == 0) ? (P0*Q+Q-1) : (rank - 1),
        rightProc = (Q0 == Q-1) ? (P0*Q) : (rank + 1);
    // blocking send and receive
#  ifdef BLOCK
    /*
    if ( rank % 2 == 1 ) {
      MPI_Send(&V(u, 1, N_loc), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm);
      MPI_Recv(&V(u, 1, 0), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, 
  	     MPI_STATUS_IGNORE);
      MPI_Send(&V(u, 1, 1), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm);
      MPI_Recv(&V(u, 1, N_loc+1), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, 
  	     comm, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&V(u, 1, 0), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, 
  	     MPI_STATUS_IGNORE);
      MPI_Send(&V(u, 1, N_loc), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm);
      MPI_Recv(&V(u, 1, N_loc+1), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, 
  	     comm, MPI_STATUS_IGNORE);
      MPI_Send(&V(u, 1, 1), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm);
    }
    */
#  endif
  // non-blocking send and receive
#  ifndef BLOCK
    double *buflef, *bufrig, *bufReclef, *bufRecrig;
    if ( w != 1 ) { 
      buflef = malloc(M_loc*w*sizeof(double)); 
      bufrig = malloc(M_loc*w*sizeof(double)); 
      bufReclef = malloc(M_loc*w*sizeof(double)); 
      bufRecrig = malloc(M_loc*w*sizeof(double)); 
      // first collect w-columns into buf, then send the buffer
      for (i = 0; i < M_loc; i++) {
        for (j = 0; j < w; j++) {
          buflef[i*w+j] = V(u,(i+w),(j+w));
  	bufrig[i*w+j] = V(u,(i+w),(N_loc+j));
        }
      }
      MPI_Isend(bufrig, M_loc*w, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqTR);
      MPI_Irecv(bufReclef, M_loc*w, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqFL);
      MPI_Isend(buflef, M_loc*w, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqTL);
      MPI_Irecv(bufRecrig, M_loc*w, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqFR);
    } else {
      MPI_Isend(&V(u, 1, N_loc), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqTR);
      MPI_Irecv(&V(u, 1, 0), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqFL);
      MPI_Isend(&V(u, 1, 1), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqTL);
      MPI_Irecv(&V(u, 1, N_loc+1), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqFR);
    }
    MPI_Wait(&reqTR, MPI_STATUS_IGNORE);
    MPI_Wait(&reqFL, MPI_STATUS_IGNORE);
    MPI_Wait(&reqTL, MPI_STATUS_IGNORE);
    MPI_Wait(&reqFR, MPI_STATUS_IGNORE);
    if ( w != 1) {
    // scatter the received buffer to halo
      for ( i = 0; i < M_loc; i++) {
        for ( j = 0; j < w; j++) {
          V(u,(w+i),j) = bufReclef[i*w+j];
          V(u,(w+i),(N_loc+w+j)) = bufRecrig[i*w+j];
        }
      }
      free(buflef);
      free(bufrig);
      free(bufReclef);
      free(bufRecrig);
    }
#  endif
  }

  //top and bottom; we can get the corner elements from the top & bottom
  if (P == 1) {
    for (j = 0; j < N_loc+2*w; j++) {
      for (i = 0; i < w; i++) {
        V(u, i, j) = V(u, M_loc+i, j);
        V(u, M_loc+w+i, j) = V(u, w+i, j);      
      }
    }
  } else {
    int topProc = (P0 == 0) ? ((P-1)*Q+Q0) : (rank - Q),
	botProc = (P0 == P-1) ? (Q0) : (rank + Q);
    double *buftop, *bufbot, *bufRectop, *bufRecbot;
    buftop = malloc((N_loc+2*w)*w*sizeof(double)); 
    bufbot = malloc((N_loc+2*w)*w*sizeof(double)); 
    bufRectop = malloc((N_loc+2*w)*w*sizeof(double)); 
    bufRecbot = malloc((N_loc+2*w)*w*sizeof(double)); 
    //first collect a row into buf, then send the buffer
    for (j = 0; j < N_loc + 2*w; j++) {
      for (i = 0; i < w; i++) {
        buftop[i*(N_loc+2*w)+j] = V(u, (w+i), j);
        bufbot[i*(N_loc+2*w)+j] = V(u, (M_loc+i), j);
      }
    }
    //communicate with top and bottom
    MPI_Isend(buftop, (N_loc+2*w)*w, MPI_DOUBLE, topProc, HALO_TAG, comm, &reqTR);
    MPI_Irecv(bufRecbot, (N_loc+2*w)*w, MPI_DOUBLE, botProc, HALO_TAG, comm, &reqFL);
    MPI_Isend(bufbot, (N_loc+2*w)*w, MPI_DOUBLE, botProc, HALO_TAG, comm, &reqTL);
    MPI_Irecv(bufRectop, (N_loc+2*w)*w, MPI_DOUBLE, topProc, HALO_TAG, comm, &reqFR);
    MPI_Wait(&reqTR, MPI_STATUS_IGNORE);
    MPI_Wait(&reqFL, MPI_STATUS_IGNORE);
    MPI_Wait(&reqTL, MPI_STATUS_IGNORE);
    MPI_Wait(&reqFR, MPI_STATUS_IGNORE);
    // scatter the received buffer to halo
    for (j = 0; j < N_loc + 2*w; j++) {
      for (i = 0; i < w; i++) {
        V(u, i, j) = bufRectop[i*(N_loc+2*w)+j];
        V(u, (M_loc+w+i), j) = bufRecbot[i*(N_loc+2*w)+j];
      }
    }
    free(buftop);
    free(bufbot);
    free(bufRectop);
    free(bufRecbot);
  }

} //updateBoundary()


// evolve advection over r timesteps, with (u,ldu) containing the local field
void parAdvect(int reps, double *u, int ldu) {
  int r, w = 1;
  for (r = 0; r < reps; r++) {
    updateBoundary(1, u, ldu);
    updateAdvectField(M_loc, N_loc, &V(u,w,w), ldu);

    if (verbosity > 2) {
      char s[64]; sprintf(s, "%d reps: u", r+1);
      printAdvectField(rank, s, M_loc+2, N_loc+2, u, ldu);
    }
  }
}

// overlap communication variant
//
// Algorithm: divide the field u into three parts.
// 1. halo part
// 2. bulk part(inner part). bulk part is the elements which don't adjacent to halo.
// 3. boundary part. it is the boundary of bulk, and it is adjacent to halo.
//
// non-blocking communication of halo -> update bulk -> MPI_Wait
// -> update boundary.
void parAdvectOverlap(int reps, double *u, int ldu) {
  int r, w = 1;
  int i, j;
  double *utmpt, *utmpb, *utmpl, *utmpr;
  int ldutmpt = 3;
  int ldutmpb = 3;
  int ldutmpl = ldu;
  int ldutmpr = ldu;
  MPI_Request reqTR, reqTL, reqFR, reqFL;
  // only work for P=1.
  if (P != 1) {
    printf("P should be 1, use -x for P != 1\n");
    exit(0);
  }
  // main loop, repetition
  for (r=0; r < reps; r++) {
    // 0. uTemp, which is temporary u field. used to calculate boundary
    utmpt = malloc(3*(N_loc+2*w)*sizeof(double)); // for top boundary
    utmpb = malloc(3*(N_loc+2*w)*sizeof(double)); // for bottom boundary
    utmpl = malloc(ldu*3*sizeof(double)); // for left boundary
    utmpr = malloc(ldu*3*sizeof(double)); // for right boundary
    for (i = 0; i < 3; i++) {
      for (j = 0; j < N_loc+2; j++) {
        V(utmpt, i, j) = V(u, i, j);
        V(utmpb, i, j) = V(u, M_loc-1+i, j);
      }
    }      
    for (i = 0; i < M_loc+2; i++) {
      for (j = 0; j < 3; j++) {
        V(utmpl, i, j) = V(u, i, j);
        V(utmpr, i, j) = V(u, i, N_loc-1+j);
      }
    }      
    // 1. update halo by nonblocking, and will not MPI_Wait at this point.
    if (Q == 1) { 
      for (i = 1; i < M_loc+1; i++) {
        V(u, i, 0) = V(u, i, N_loc);
        V(u, i, N_loc+1) = V(u, i, 1);
      }
    } else {
      int leftProc = (Q0 == 0) ? (P0*Q+Q-1) : (rank - 1),
          rightProc = (Q0 == Q-1) ? (P0*Q) : (rank + 1);
      MPI_Isend(&V(u, 1, N_loc), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqTR);
      MPI_Irecv(&V(u, 1, 0), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqFL);
      MPI_Isend(&V(u, 1, 1), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqTL);
      MPI_Irecv(&V(u, 1, N_loc+1), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqFR);
    }
    // 2. calculate bulk
    updateAdvectField(M_loc-2, N_loc-2,&V(u,(w+1),(w+1)), ldu); //ldu never used
    // 3. calculate boundary
    //  3.1 wait
    if (Q != 1) {
      MPI_Wait(&reqTR, MPI_STATUS_IGNORE);
      MPI_Wait(&reqFL, MPI_STATUS_IGNORE);
      MPI_Wait(&reqTL, MPI_STATUS_IGNORE);
      MPI_Wait(&reqFR, MPI_STATUS_IGNORE);
    }
    //  3.2 set top bottom of u, and the halo part of utmpt, utmpb, utmpr, utmpl
    for ( i =  1; i < M_loc+1; i++) {
      V(utmpl, i, 0) = V(u, i, 0);
      V(utmpr, i, 2) = V(u, i, N_loc+1);
    }
    for ( i =  1; i < 3; i++) {
      V(utmpt, i, 0) = V(u, i, 0);
      V(utmpt, i, N_loc+1) = V(u, i, N_loc+1);
    }
    for ( i =  0; i < 2; i++) {
      V(utmpb, i, 0) = V(u, M_loc-1+i, 0);
      V(utmpb, i, N_loc+1) = V(u, M_loc-1+i, N_loc+1);
    }
    if (P == 1) {
      for (j = 0; j < N_loc+2; j++) {
        V(u, 0, j) = V(u, M_loc, j);
        V(u, M_loc+1, j) = V(u, 1, j);      
        V(utmpt, 0, j) = V(utmpb, 1, j);
        V(utmpb, 2, j) = V(utmpt, 1, j);      
      }
      for (j = 0; j < 3; j++) { 
        V(utmpl, 0, j) = V(utmpl, M_loc, j);
        V(utmpl, M_loc+1, j) = V(utmpl, 1, j);      
        V(utmpr, 0, j) = V(utmpr, M_loc, j);
        V(utmpr, M_loc+1, j) = V(utmpr, 1, j);      
      }
    } 
    //  3.3 calculate top boundary. using utmpt !!!
    updateAdvectField(1,N_loc, &V(utmpt,1,1), ldutmpt);
    //  3.4 calculate left boundary. using utmpl !!!
    updateAdvectField(M_loc-2,1, &V(utmpl,2,1), ldutmpl); 
    //  3.5 calculate right boundary. using utmpr !!!
    updateAdvectField(M_loc-2,1, &V(utmpr,2,1), ldutmpr);
    //  3.6 calculate bottom boundary. using utmpb !!!
    updateAdvectField(1,N_loc, &V(utmpb,1,1), ldutmpb);
    // 4. copy boundary of uTemp to u
    for ( i =  1; i < M_loc+1; i++) { // left and right boundary
      V(u,i,1) = V(utmpl,i,1);
      V(u,i,N_loc) = V(utmpr,i,1);
    }
    for ( j = 1; j < N_loc+1; j++) { // top and bottom boundary
      V(u,1,j) = V(utmpt,1,j);
      V(u,M_loc,j) = V(utmpb,1,j);
    }
    free(utmpb);
    free(utmpt);
    free(utmpl);
    free(utmpr);
  } // r loop
} //parAdvectOverlap

// wide halo variant
void parAdvectWide(int reps, int w, double *u, int ldu) {
  int tmp, tmp2, r;
  for (r = 0; r < reps; r++) {
    if (r % w == 0) {
      updateBoundary(w, u, ldu);
    }
    tmp = r % w;
    tmp2 = tmp + 1;
    updateAdvectField(M_loc+2*(w-1-tmp), N_loc+2*(w-1-tmp), &V(u,tmp2,tmp2), ldu);
  }
}

// extra optimization variant
void parAdvectExtra(int reps, double *u, int ldu) {
  int r, w = 1;
  int i, j;
  double *utmpt, *utmpb, *utmpl, *utmpr;
  int ldutmpt = 3;
  int ldutmpb = 3;
  int ldutmpl = ldu;
  int ldutmpr = ldu;
  MPI_Request reqTR, reqTL, reqFR, reqFL;
  // main loop, repetition
  for (r=0; r < reps; r++) {
    // 0. uTemp, which is temporary u field. used to calculate boundary
    utmpt = malloc(3*(N_loc+2*w)*sizeof(double)); // for top boundary
    utmpb = malloc(3*(N_loc+2*w)*sizeof(double)); // for bottom boundary
    utmpl = malloc(ldu*3*sizeof(double)); // for left boundary
    utmpr = malloc(ldu*3*sizeof(double)); // for right boundary
    for (i = 0; i < 3; i++) {
      for (j = 0; j < N_loc+2; j++) {
        V(utmpt, i, j) = V(u, i, j);
        V(utmpb, i, j) = V(u, M_loc-1+i, j);
      }
    }      
    for (i = 0; i < M_loc+2; i++) {
      for (j = 0; j < 3; j++) {
        V(utmpl, i, j) = V(u, i, j);
        V(utmpr, i, j) = V(u, i, N_loc-1+j);
      }
    }      
    // 1. update left and right halo by nonblocking, and will not MPI_Wait at this point.
    if (Q == 1) { 
      for (i = 1; i < M_loc+1; i++) {
        V(u, i, 0) = V(u, i, N_loc);
        V(u, i, N_loc+1) = V(u, i, 1);
      }
    } else {
      int leftProc = (Q0 == 0) ? (P0*Q+Q-1) : (rank - 1),
          rightProc = (Q0 == Q-1) ? (P0*Q) : (rank + 1);
      MPI_Isend(&V(u, 1, N_loc), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqTR);
      MPI_Irecv(&V(u, 1, 0), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqFL);
      MPI_Isend(&V(u, 1, 1), M_loc, MPI_DOUBLE, leftProc, HALO_TAG, comm, &reqTL);
      MPI_Irecv(&V(u, 1, N_loc+1), M_loc, MPI_DOUBLE, rightProc, HALO_TAG, comm, &reqFR);
    }
    // 2. calculate bulk
    updateAdvectField(M_loc-2, N_loc-2,&V(u,(w+1),(w+1)), ldu); //ldu never used
    // 3. wait for left right halo
    if (Q != 1) {
      MPI_Wait(&reqTR, MPI_STATUS_IGNORE);
      MPI_Wait(&reqFL, MPI_STATUS_IGNORE);
      MPI_Wait(&reqTL, MPI_STATUS_IGNORE);
      MPI_Wait(&reqFR, MPI_STATUS_IGNORE);
    }
    // 4. update the relevant element in utmpb, utmpt, utmpl, utmpr
    for ( i =  1; i < M_loc+1; i++) {
      V(utmpl, i, 0) = V(u, i, 0);
      V(utmpr, i, 2) = V(u, i, N_loc+1);
    }
    for ( i =  1; i < 3; i++) {
      V(utmpt, i, 0) = V(u, i, 0);
      V(utmpt, i, N_loc+1) = V(u, i, N_loc+1);
    }
    for ( i =  0; i < 2; i++) {
      V(utmpb, i, 0) = V(u, M_loc-1+i, 0);
      V(utmpb, i, N_loc+1) = V(u, M_loc-1+i, N_loc+1);
    }
    // 5. update top and bottom halo by nonblocking, and do MPI_WAIT immediatly 
    if (P == 1) {
      for (j = 0; j < N_loc+2; j++) {
        V(u, 0, j) = V(u, M_loc, j);
        V(u, M_loc+1, j) = V(u, 1, j);      
        V(utmpt, 0, j) = V(utmpb, 1, j);
        V(utmpb, 2, j) = V(utmpt, 1, j);      
      }
      for (j = 0; j < 3; j++) { 
        V(utmpl, 0, j) = V(utmpl, M_loc, j);
        V(utmpl, M_loc+1, j) = V(utmpl, 1, j);      
        V(utmpr, 0, j) = V(utmpr, M_loc, j);
        V(utmpr, M_loc+1, j) = V(utmpr, 1, j);      
      }
    } else {
      int topProc = (P0 == 0) ? ((P-1)*Q+Q0) : (rank - Q),
  	botProc = (P0 == P-1) ? (Q0) : (rank + Q);
      double *buftop, *bufbot, *bufRectop, *bufRecbot;
      buftop = malloc((N_loc+2*w)*w*sizeof(double)); 
      bufbot = malloc((N_loc+2*w)*w*sizeof(double)); 
      bufRectop = malloc((N_loc+2*w)*w*sizeof(double)); 
      bufRecbot = malloc((N_loc+2*w)*w*sizeof(double)); 
      //first collect a row into buf, then send the buffer
      for (j = 0; j < N_loc + 2*w; j++) {
          buftop[j] = V(utmpt, 1, j);
          bufbot[j] = V(utmpb, 1, j);
      }
      //communicate with top and bottom
      MPI_Isend(buftop, (N_loc+2*w)*w, MPI_DOUBLE, topProc, HALO_TAG, comm, &reqTR);
      MPI_Irecv(bufRecbot, (N_loc+2*w)*w, MPI_DOUBLE, botProc, HALO_TAG, comm, &reqFL);
      MPI_Isend(bufbot, (N_loc+2*w)*w, MPI_DOUBLE, botProc, HALO_TAG, comm, &reqTL);
      MPI_Irecv(bufRectop, (N_loc+2*w)*w, MPI_DOUBLE, topProc, HALO_TAG, comm, &reqFR);
      MPI_Wait(&reqTR, MPI_STATUS_IGNORE);
      MPI_Wait(&reqFL, MPI_STATUS_IGNORE);
      MPI_Wait(&reqTL, MPI_STATUS_IGNORE);
      MPI_Wait(&reqFR, MPI_STATUS_IGNORE);
      // scatter the received buffer to halo
      for (j = 0; j < N_loc + 2*w; j++) {
          V(u, 0, j) = bufRectop[j];
          V(u, M_loc+1, j) = bufRecbot[j];
          V(utmpt, 0, j) = bufRectop[j];
          V(utmpb, 2, j) = bufRecbot[j];
      }
      for ( j =  0; j < 3; j++) {
        V(utmpl, 0, j) = V(u, 0, j);
        V(utmpl, M_loc+1, j) = V(u, M_loc+1, j);
        V(utmpr, 0, j) = V(u, 0, N_loc-1+j);
        V(utmpr, M_loc+1, j) = V(u, M_loc+1, N_loc-1+j);
      }
      free(buftop);
      free(bufbot);
      free(bufRectop);
      free(bufRecbot);
    }
    // 6. calculate boundary
    //  calculate top boundary. using utmpt !!!
    updateAdvectField(1,N_loc, &V(utmpt,1,1), ldutmpt);
    //  calculate left boundary. using utmpl !!!
    updateAdvectField(M_loc-2,1, &V(utmpl,2,1), ldutmpl); 
    //  calculate right boundary. using utmpr !!!
    updateAdvectField(M_loc-2,1, &V(utmpr,2,1), ldutmpr);
    //  calculate bottom boundary. using utmpb !!!
    updateAdvectField(1,N_loc, &V(utmpb,1,1), ldutmpb);
    // 7. copy boundary of uTemp to u
    for ( i =  1; i < M_loc+1; i++) { // left and right boundary
      V(u,i,1) = V(utmpl,i,1);
      V(u,i,N_loc) = V(utmpr,i,1);
    }
    for ( j = 1; j < N_loc+1; j++) { // top and bottom boundary
      V(u,1,j) = V(utmpt,1,j);
      V(u,M_loc,j) = V(utmpb,1,j);
    }
    free(utmpb);
    free(utmpt);
    free(utmpl);
    free(utmpr);
  } // r loop
}
