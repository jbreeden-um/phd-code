#include "verify.h"

/* Generate all the test quaternions at once */
/* These two functions work by projecting an even grid of a unit hypercube onto a unit hypersphere */
int generate_unit_quaternions(double sparsity, double **Q){
    int n_per_dimension = (int)(2.0/sparsity+0.5) + 1;
    double increment = 2.0/(n_per_dimension-1);
    int length = 8*n_per_dimension*n_per_dimension*n_per_dimension;
    int count = 0;

    // printf("Expected: %d\n", length);
    *Q = malloc(4*length*sizeof(double));
    
    double w, x, y, z, mag;
    int d0, d1, d2, d3;
    for (int fixed_dimension=0; fixed_dimension<4; fixed_dimension++){
        d0 = fixed_dimension;
        d1 = (fixed_dimension + 1) % 4;
        d2 = (fixed_dimension + 2) % 4;
        d3 = (fixed_dimension + 3) % 4;
        for (int side=-1; side<=1; side+=2){
            w = (double)side;
            for (int i=0; i<n_per_dimension; i++){
                x = -1.0 + increment*i;
                for (int j=0; j<n_per_dimension; j++){
                    y = -1.0 + increment*j;
                    for (int k=0; k<n_per_dimension; k++){
                        z = -1.0 + increment*k;
                        mag = sqrt(w*w + x*x + y*y + z*z);
                        (*Q)[4*count+d0] = w/mag;
                        (*Q)[4*count+d1] = x/mag;
                        (*Q)[4*count+d2] = y/mag;
                        (*Q)[4*count+d3] = z/mag;
                        count++;
                    }
                }
            }
        }
    }
    // printf("Points: %d\n", count);
    // printf("Matches: %d\n", count == length);
    return length;
}

/* Generate one quaternion at a time to save memory */
int generate_next_quaternion(double sparsity, double Q[4], char first){
    static int n_per_dimension = 2;
    static double increment = 0.1;
    
    static int curr_fixed_dimension = 0;
    static int curr_side = -1;
    static int curr_i = 0;
    static int curr_j = 0;
    static int curr_k = 0;

    if (first) {
        n_per_dimension = (int)(2.0/sparsity+0.5) + 1;
        increment = 2.0/(n_per_dimension-1);
        curr_fixed_dimension = 0;
        curr_side = -1;
        curr_i = 0;
        curr_j = 0;
        curr_k = 0;
        int length = 7*n_per_dimension*n_per_dimension*n_per_dimension;
        return length;
    }
    
    double w, x, y, z, mag;
    int d0, d1, d2, d3;
    for (int fixed_dimension=curr_fixed_dimension; fixed_dimension<4; fixed_dimension++){
        d0 = fixed_dimension;
        d1 = (fixed_dimension + 1) % 4;
        d2 = (fixed_dimension + 2) % 4;
        d3 = (fixed_dimension + 3) % 4;
        for (int side=curr_side; side<=1; side+=2){
            w = (double)side;
            for (int i=curr_i; i<n_per_dimension; i++){
                x = -1.0 + increment*i;
                for (int j=curr_j; j<n_per_dimension; j++){
                    y = -1.0 + increment*j;
                    for (int k=curr_k; k<n_per_dimension; k++){
                        z = -1.0 + increment*k;
                        mag = sqrt(w*w + x*x + y*y + z*z);
                        Q[d0] = w/mag;
                        Q[d1] = x/mag;
                        Q[d2] = y/mag;
                        Q[d3] = z/mag;
                        curr_k++;
                        return 1;
                    }
                    curr_j++;
                    curr_k = 0;
                }
                curr_i++;
                curr_j = 0;
            }
            curr_side+=2;
            curr_i = 0;
        }
        curr_fixed_dimension++;
        curr_side = -1;
        if (curr_fixed_dimension == 3) curr_side = 1; // exclude negative rotation quaternions
    }
    
    printf("Something is wrong in generate_next_quaternion.\n");
    return -1;
}

/* Given a quaternion and one CBF, generate several w values that result in the CBF being equal to zero */
int generate_w1(double q[4], double M[4][4], double sparsity, int Nmax, double *out){
    double h = h_func(q, M);
    if (h > 0){
        printf("Invalid point. h is positive");
        exit(1);
    }
    double b = sqrt(-2*mu*h);
    double A[3];
    hdot_affine(q, M, A);

    double n = Mag3(A);
    if (n < eps_tol){
        printf("hdot sensitivity matrix is zero in solve_for_w1. You should look into that.\n");
    }

    double par[3];
    for (int i=0; i<3; i++) par[i] = A[i]*b/(n*n);
    
    // construct basis for the w nullspace
    double orth[2][3];
    double ei[3] = {0.0,0.0,0.0};
    int pivot;
    if (fabs(A[0]) < fabs(A[1])) pivot = 0;
    else pivot = 1;
    if (fabs(A[2]) < fabs(A[pivot])) pivot = 2;
    ei[pivot] = 1.0;
    Cross3(par, ei, orth[0]);
    n = Mag3(orth[0]);
    for (int i=0; i<3; i++) orth[0][i] /= n;
    Cross3(par, orth[0], orth[1]);
    n = Mag3(orth[1]);
    for (int i=0; i<3; i++) orth[1][i] /= n;

    int n_per_dimension = MIN( (int)(2.0*w_max/sparsity+0.5) + 1, (int)sqrt(Nmax/sqrt(2.0)) );
    double increment = 2.0*w_max/(n_per_dimension-1);
    int count = 0;
    
    double vec[3], x, y;
    if (Mag3(par) > w_max*sqrt(3)){
        return 0;
    }
    for (int i=0; i<2*n_per_dimension-1; i++){
        // Start from -2.0*w_max instead of -w_max to ensure all valid points are covered when orth is added to par
        x = -2.0*w_max + increment*i;
        for (int j=0; j<2*n_per_dimension-1; j++){
            y = -2.0*w_max + increment*j;
            for (int p=0; p<3; p++) vec[p] = par[p] + x*orth[0][p] + y*orth[1][p];
            if (fabs(vec[0]) <= w_max && fabs(vec[1]) <= w_max && fabs(vec[2]) <= w_max){
                if (count > Nmax){
                    printf("Number of w points exceeds maximum allocated memory in generate_w1. You should fix that. Exiting\n");
                    exit(1);   
                }
                for (int p=0; p<3; p++) out[3*count+p] = vec[p];
                count++;
            }
        }
    }

    return count;
}

/* Given a quaternion and two CBFs, solve for the space of w values that result in both CBFs being zero, and return a sampling of that space */
int solve_for_w2(double q[4], double M1[4][4], double M2[4][4], int Nmin, int Nmax, double sparsity, double *out){
    double h1 = h_func(q, M1);
    double h2 = h_func(q, M2);
    if (h1 > 0 || h2 > 0){
        printf("Invalid point. h is positive");
        exit(1);
    }
    double b[2] = {sqrt(-2*mu*h1), sqrt(-2*mu*h2)};
    double A1[3], A2[3];
    hdot_affine(q, M1, A1);
    hdot_affine(q, M2, A2);

    if (Mag3(A1) < eps_tol || Mag3(A2) < eps_tol){
        if (h1 > -eps_tol && h2 > -eps_tol) 
            printf("hdot sensitivity matrix is zero in solve_for_w2. You should look into that.\n");
        else
            // for this system, points directly opposite the constraint have zero sensitivity.
            // thus, if h != 0, there is no value of w that will cause H to equal zero.
            return 0;
    }
    
    double orth[3];
    Cross3(A1, A2, orth);

    double par[3];
    if (Mag3(orth)/(Mag3(A1)*Mag3(A2)) < eps_tol){
        if (fabs(Mag3(A1)/Mag3(A2) - fabs(b[0]/b[1])) < eps_tol && Dot3(A1, A2)*b[0]*b[1] >= 0){
            printf("The pseudoinverse is degenerate, but a solution still exists. You should look into that. Returning empty for now.\n");
            printf("The current point is q = [%lf  %lf  %lf  %lf]\n", q[0], q[1], q[2], q[3]);
        } else {
            // The pseudo-inverse is degenerate, but it's okay because no solution exists anyways.
            return 0;
        }
    }
    pseudo_inverse_2x3(A1, A2, b, par);
    if (isnan(par[0]) || isnan(par[1]) || isnan(par[2])){
        printf("The pseudoinverse is returning NaNs. You should look into that.\n");
    }

    double alpha1, alpha2, alpha;
    double alpha_allow_min[3], alpha_allow_max[3];
    for (int i=0; i<3; i++){
        if (fabs(orth[i]) > eps_tol){
            alpha1 = (w_max - par[i])/orth[i];
            alpha2 = (-w_max - par[i])/orth[i];
            if (par[i] <= w_max && par[i] >= -w_max){
                alpha_allow_max[i] = MAX(alpha1, alpha2);
                alpha_allow_min[i] = MIN(alpha1, alpha2);
            } else if (par[i] > 0){ // par[i] is too large and needs to decrease
                if (orth[i] < 0){ // positive alpha1 brings par back down to the acceptable range
                    alpha_allow_min[i] = alpha1; // alpha must be at least this large
                    alpha_allow_max[i] = alpha2;
                } else { // negative alpha1 brings par back down to the acceptable range
                    alpha_allow_max[i] = alpha1; // alpha must be at most this large
                    alpha_allow_min[i] = alpha2;
                }
            } else { // par[i] is too small and needs to increase
                if (orth[i] > 0){ // positive alpha2 brings par back up to the acceptable range
                    alpha_allow_min[i] = alpha2; // alpha must be at least this large
                    alpha_allow_max[i] = alpha1;
                } else { // negative alpha2 brings par back up to the acceptable range
                    alpha_allow_max[i] = alpha2; // alpha must be at most this large
                    alpha_allow_min[i] = alpha1;
                }
            }
        } else {
            if (par[i] <= w_max && par[i] >= -w_max){
                alpha_allow_max[i] = 1e10;
                alpha_allow_min[i] = -1e10;
            } else {
                // the required angular velocity is unachievable and no addition from orth will change that
                alpha_allow_max[i] = -1e10;
                alpha_allow_min[i] = 1e10;
            }
        }
    }
    double alpha_min = MAX(alpha_allow_min[0], MAX(alpha_allow_min[1], alpha_allow_min[2]));
    double alpha_max = MIN(alpha_allow_max[0], MIN(alpha_allow_max[1], alpha_allow_max[2]));
    if (alpha_min > alpha_max) return 0;

    if (Nmin > Nmax || Nmin < 1 || Nmax < 1){
        printf("Invalid requested dimensions.\n");
        exit(1);
    }

    int Ndes = (int)(2.0*w_max/sparsity)+1;
    int N = Ndes;
    if (N > Nmax) N = Nmax;
    if (N < Nmin) N = Nmin;
    double increment = (alpha_max - alpha_min)/(N-1);
    if (N==1) increment = 0.0;
    for (int i=0; i<N; i++){
        alpha = alpha_min + i*increment;
        for (int j=0; j<3; j++){
            out[3*i+j] = par[j] + orth[j]*alpha;
            if (isnan(out[3*i+j])){
                printf("Huh? There is a nan in the output of solve_for_w2.\n");
            }
        }
    }
    return N;
}