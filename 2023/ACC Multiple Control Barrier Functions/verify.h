#ifndef __VERIFY__
#define __VERIFY__

#ifdef _RUN_CASE_1_
    #undef _RUN_CASE_1_
#endif
#ifdef _RUN_CASE_2_
    #undef _RUN_CASE_2_
#endif
#ifdef _RUN_CASE_3_
    #undef _RUN_CASE_3_
#endif

// #define _RUN_CASE_1_   // this is a sparse run
#define _RUN_CASE_2_   // this is the main run
// #define _RUN_CASE_3_   // this is a run without clustering

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define NM0 (2)
#define MAX_CONSTRAINTS (12)
#define INIT_TRACKED_POINTS (100000)
#define K_MAX_CLUSTER (20)
extern long MAX_TRACKED_POINTS;
extern const double q_neighbor_tol;

extern double q_sparsity;
extern double w_sparsity;
extern const int Nwmin;
extern const int Nwmax;

extern const double mu;
extern const double w_max;
extern const double u_max;
extern const double H_zero_tol;
extern const double eps_tol;
extern double rb[3];
extern double M_start[NM0][4][4];
extern double (*M_all)[4][4];
extern int NM;

double h_func(double q[4], double M[4][4]);
double hdot_func(double q[4], double w[3], double M[4][4]);
double H_func(double q[4], double w[3], double M[4][4]);
void hdot_affine(double q[4], double M[4][4], double out[3]);
double hddot_func(double q[4], double w[3], double u[3], double M[4][4]);
double Hdot_func(double q[4], double w[3], double u[3], double M[4][4]);

double Det3(double A[3][3]);
void Inv3(double A[3][3], double b[3], double x[3]);
void pseudo_inverse_2x3(double row1[3], double row2[3], double b[2], double out[3]);
double MIN(double a, double b);
double MAX(double a, double b);
double Dot3(double a[3], double b[3]);
double Mag3(double a[3]);
void Cross3(double A1[3], double A2[3], double orth[3]);
void RotQ(double u[3], double q[4], double out[3]);
void QxQ(double q1[4], double q2[3], double out[4]);
void QxQT(double q1[4], double q2[3], double out[4]);

int generate_unit_quaternions(double sparsity, double **Q);
int generate_next_quaternion(double sparsity, double Q[4], char first);
int generate_w1(double q[4], double M[4][4], double sparsity, int Nmax, double *out);
int solve_for_w2(double q[4], double M1[4][4], double M2[4][4], int Nmin, int Nmax, double sparsity, double *out);

int check_points(double (**track)[8], char quit);

#endif