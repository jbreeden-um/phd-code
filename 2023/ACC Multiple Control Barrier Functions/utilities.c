#include "verify.h"

double Det3(double A[3][3]){
    return A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0];
}

void Inv3(double A[3][3], double b[3], double x[3]){
    double d = Det3(A);
    x[0] =  (A[0][1]*A[1][2]*b[2] - A[0][2]*A[1][1]*b[2] - A[0][1]*A[2][2]*b[1] + A[0][2]*A[2][1]*b[1] + A[1][1]*A[2][2]*b[0] - A[1][2]*A[2][1]*b[0])/d;
    x[1] = -(A[0][0]*A[1][2]*b[2] - A[0][2]*A[1][0]*b[2] - A[0][0]*A[2][2]*b[1] + A[0][2]*A[2][0]*b[1] + A[1][0]*A[2][2]*b[0] - A[1][2]*A[2][0]*b[0])/d;
    x[2] =  (A[0][0]*A[1][1]*b[2] - A[0][1]*A[1][0]*b[2] - A[0][0]*A[2][1]*b[1] + A[0][1]*A[2][0]*b[1] + A[1][0]*A[2][1]*b[0] - A[1][1]*A[2][0]*b[0])/d;
}

void pseudo_inverse_2x3(double row1[3], double row2[3], double b[2], double out[3]){
    double x1 = row1[0];
    double x2 = row1[1];
    double x3 = row1[2];
    double y1 = row2[0];
    double y2 = row2[1];
    double y3 = row2[2];
    double z1 = b[0];
    double z2 = b[1];
    out[0] = - z2*((x1*(x1*y1 + x2*y2 + x3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2) - (y1*(x1*x1 + x2*x2 + x3*x3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2)) - z1*((y1*(x1*y1 + x2*y2 + x3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2) - (x1*(y1*y1 + y2*y2 + y3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2));
    out[1] = - z2*((x2*(x1*y1 + x2*y2 + x3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2) - (y2*(x1*x1 + x2*x2 + x3*x3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2)) - z1*((y2*(x1*y1 + x2*y2 + x3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2) - (x2*(y1*y1 + y2*y2 + y3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2));
    out[2] = - z2*((x3*(x1*y1 + x2*y2 + x3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2) - (y3*(x1*x1 + x2*x2 + x3*x3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2)) - z1*((y3*(x1*y1 + x2*y2 + x3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2) - (x3*(y1*y1 + y2*y2 + y3*y3))/(x1*x1*y2*y2 + x1*x1*y3*y3 - 2*x1*x2*y1*y2 - 2*x1*x3*y1*y3 + x2*x2*y1*y1 + x2*x2*y3*y3 - 2*x2*x3*y2*y3 + x3*x3*y1*y1 + x3*x3*y2*y2));
}

double MIN(double a, double b){
    if (a < b) return a;
    else return b;
}

double MAX(double a, double b){
    if (a > b) return a;
    else return b;
}

double Dot3(double a[3], double b[3]){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double Mag3(double a[3]){
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

void Cross3(double A1[3], double A2[3], double orth[3]){
    orth[0] = A1[1]*A2[2] - A1[2]*A2[1];
    orth[1] = A1[2]*A2[0] - A1[0]*A2[2];
    orth[2] = A1[0]*A2[1] - A1[1]*A2[0];
}

void RotQ(double u[3], double q[4], double out[3]){
    out[0] = q[0]*(q[0]*u[0] + q[1]*u[1] + q[2]*u[2]) + q[1]*(q[0]*u[1] - q[1]*u[0] + q[3]*u[2]) - q[2]*(q[2]*u[0] - q[0]*u[2] + q[3]*u[1]) + q[3]*(q[1]*u[2] - q[2]*u[1] + q[3]*u[0]);
    out[1] = q[1]*(q[0]*u[0] + q[1]*u[1] + q[2]*u[2]) - q[0]*(q[0]*u[1] - q[1]*u[0] + q[3]*u[2]) + q[2]*(q[1]*u[2] - q[2]*u[1] + q[3]*u[0]) + q[3]*(q[2]*u[0] - q[0]*u[2] + q[3]*u[1]);
    out[2] = q[2]*(q[0]*u[0] + q[1]*u[1] + q[2]*u[2]) + q[0]*(q[2]*u[0] - q[0]*u[2] + q[3]*u[1]) - q[1]*(q[1]*u[2] - q[2]*u[1] + q[3]*u[0]) + q[3]*(q[0]*u[1] - q[1]*u[0] + q[3]*u[2]);
}

void QxQ(double q1[4], double q2[3], double out[4]){
    out[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
    out[1] = q1[3]*q2[1] - q1[0]*q2[2] + q1[1]*q2[3] + q1[2]*q2[0];
    out[2] = q1[3]*q2[2] + q1[0]*q2[1] - q1[1]*q2[0] + q1[2]*q2[3];
    out[3] = q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2];
}

void QxQT(double q1[4], double q2[3], double out[4]){
    out[0] = -q1[3]*q2[0] + q1[0]*q2[3] - q1[1]*q2[2] + q1[2]*q2[1];
    out[1] = -q1[3]*q2[1] + q1[0]*q2[2] + q1[1]*q2[3] - q1[2]*q2[0];
    out[2] = -q1[3]*q2[2] - q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3];
    out[3] =  q1[3]*q2[3] + q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2];
}
