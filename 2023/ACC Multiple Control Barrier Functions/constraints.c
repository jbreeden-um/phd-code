#include "verify.h"

double h_func(double q[4], double M[4][4]){
    double out = 0.0;
    for (int i=0; i<4; i++){
        for (int j=0; j<4;  j++){
            out += q[i]*M[i][j]*q[j];
        }
    }
    return out;
}

double hdot_func(double q[4], double w[3], double M[4][4]){
    double W[4][4] = {{0, w[2], -w[1], w[0]},
                      {-w[2], 0, w[0], w[1]},
                      {w[1], -w[0], 0, w[2]},
                      {-w[0], -w[1], -w[2], 0}};
    
    double mat[4][4] = {0};
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            for (int k=0; k<4; k++){
                mat[i][j] += W[k][i]*M[k][j] + M[i][k]*W[k][j];
            }
        }
    }

    double out = 0.0;
    for (int i=0; i<4; i++){
        for (int j=0; j<4;  j++){
            out += q[i]*mat[i][j]*q[j]/2;
        }
    }
    return out;
}

double hddot_func(double q[4], double w[3], double u[3], double M[4][4]){
    double W[4][4] = {{0, w[2], -w[1], w[0]},
                      {-w[2], 0, w[0], w[1]},
                      {w[1], -w[0], 0, w[2]},
                      {-w[0], -w[1], -w[2], 0}};
    double U[4][4] = {{0, u[2], -u[1], u[0]},
                      {-u[2], 0, u[0], u[1]},
                      {u[1], -u[0], 0, u[2]},
                      {-u[0], -u[1], -u[2], 0}};
    
    double mat_i[4][4] = {0};
    double mat_u[4][4] = {0};
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            for (int k=0; k<4; k++){
                mat_i[i][j] += W[k][i]*M[k][j] + M[i][k]*W[k][j];
                mat_u[i][j] += U[k][i]*M[k][j] + M[i][k]*U[k][j];
            }
        }
    }

    double mat_o[4][4] = {0};
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            for (int k=0; k<4; k++){
                mat_o[i][j] += W[k][i]*mat_i[k][j] + mat_i[i][k]*W[k][j];
            }
        }
    }

    double out = 0.0;
    for (int i=0; i<4; i++){
        for (int j=0; j<4;  j++){
            out += q[i]*mat_o[i][j]*q[j]/4 + q[i]*mat_u[i][j]*q[j]/2;
        }
    }
    return out;
}

double H_func(double q[4], double w[3], double M[4][4]){
    double h = h_func(q, M);
    if (h <= 0)
        return hdot_func(q, w, M) - sqrt(-2*mu*h);
    else{
        printf("Positive h encountered. Returning +1\n");
        return 1.0;
    }
}

double Hdot_func(double q[4], double w[3], double u[3], double M[4][4]){
    double h = h_func(q, M);
    if (h < 0){
        return hddot_func(q, w, u, M) + sqrt(mu)*hdot_func(q, w, M)/sqrt(-2*h);
    } else if (h==0){
        if (fabs(H_func(q, w, M)) < H_zero_tol) return hddot_func(q, w, u, M);
        else {
            printf("Hdot is undefined in this case. You should fix that. Returning +1\n");
            return 1.0;
        }
    } else{
        printf("Positive h encountered. Returning +1\n");
        return 1.0;
    }
}

void hdot_affine(double q[4], double M[4][4], double out[3]){
    double Q[4][3] = {{q[3], -q[2], q[1]},
                      {q[2], q[3], -q[0]},
                      {-q[1], q[0], q[3]},
                      {-q[0], -q[1], -q[2]}};
    
    double mat[4][3] = {0};
    for (int i=0; i<4; i++){
        for (int j=0; j<3; j++){
            for (int k=0; k<4; k++){
                mat[i][j] += M[k][i]*Q[k][j] + M[i][k]*Q[k][j];
            }
        }
    }

    for (int i=0; i<3; i++){
        out[i] = 0.0;
        for (int j=0; j<4; j++){
            out[i] += q[j]*mat[j][i]/2;
        }
    }
}