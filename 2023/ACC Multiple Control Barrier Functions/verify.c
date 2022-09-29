#include "verify.h"

void A_func(double q[4], double M[4][4], double out[3]){
    hdot_affine(q, M, out);
}

/* b func when the CBF is equal to zero */
double b0_func(double q[4], double w[3], double M[4][4]){
    double h = h_func(q, M);
    double hdot = hdot_func(q, w, M);
    if (h >= 0){
        printf("Positive h encountered. Quitting\n");
        exit(1);
    }

    double W[4][4] = {{0, w[2], -w[1], w[0]},
                      {-w[2], 0, w[0], w[1]},
                      {w[1], -w[0], 0, w[2]},
                      {-w[0], -w[1], -w[2], 0}};
    
    double mat_i[4][4] = {0};
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            for (int k=0; k<4; k++){
                mat_i[i][j] += W[k][i]*M[k][j] + M[i][k]*W[k][j];
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
        for (int j=0; j<4; j++){
            out += q[i]*mat_o[i][j]*q[j]/4;
        }
    }
    if (h != 0)
        out += sqrt(mu)*hdot/sqrt(-2*h);

    return -out;
}

int get_active_constraints(double q[4], double w[3], double *A, double *b){
    int count = 0;
    double H;
    for (int iH=0; iH<NM; iH++){
        H = H_func(q, w, M_all[iH]);
        if (fabs(H) < H_zero_tol){
            if (count >= MAX_CONSTRAINTS){
                printf("Out of memory. Please change MAX_CONSTRAINTS");
                exit(1);
            }
            A_func(q, M_all[iH], &A[3*count]);
            b[count] = b0_func(q, w, M_all[iH]);
            count++;
        } else if (H > H_zero_tol){
            return -1; // we are outside the safe set, so do not bother going further
        }
    }

    if (count < 2) printf("Something is wrong. There should always be exactly two constraints in this case.\n");
    return count;
}


char is_viable_exact(double *A, double *b, int N)  {
    double LHS[MAX_CONSTRAINTS+6][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1}};
    double RHS[MAX_CONSTRAINTS+6] = {u_max, u_max, u_max, u_max, u_max, u_max};

    for (int i=0; i<N; i++){
        for (int j=0; j<3; j++){
            LHS[i+6][j] = A[3*i+j];
        }
        RHS[i+6] = b[i];
    }

    double mat[3][3];
    double y[3];
    double u[3];
    for (int i=0; i<N+6; i++){
        for (int j=i+1; j<N+6; j++){
            for (int k=j+1; k<N+6; k++){
                memcpy(&mat[0][0], &LHS[i][0], 3*sizeof(double));
                memcpy(&mat[1][0], &LHS[j][0], 3*sizeof(double));
                memcpy(&mat[2][0], &LHS[k][0], 3*sizeof(double));
                y[0] = RHS[i];
                y[1] = RHS[j];
                y[2] = RHS[k];
                if (fabs(Det3(mat)) > eps_tol){
                    Inv3(mat, y, u);
                    char meets_all_constraints = 1;
                    for (int m=0; m<N+6; m++){
                        if (LHS[m][0]*u[0] + LHS[m][1]*u[1] + LHS[m][2]*u[2] - RHS[m] > eps_tol){
                            meets_all_constraints = 0;
                            break;
                        }
                    }
                    if (meets_all_constraints){
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

int check_points(double (**ptr)[8], char quit){
    static int Nq = -1;
    double Q[4];
    static double *W;
    if (Nq == -1){
        W = malloc(3*Nwmax*sizeof(double));
    }

    double (*track)[8] = *ptr;

    int Nw, Ncon;
    double A[3*MAX_CONSTRAINTS];
    double b[MAX_CONSTRAINTS];
    double *h = malloc(NM*sizeof(double));
    double *q, *w;
    char inside_S;
    char valid_point;
    int n_points = 0;
    int n_opt = 0;
    int n_qvalid = 0;

    int N_track = 0;

    char output[20];
    sprintf(output, "data/conflicts%d.dat", NM);
    FILE *fout = fopen(output, "w+");
    int combinations_tried = 0;
    for (int c1=0; c1<NM; c1++){
        for (int c2=c1+1; c2<NM; c2++){
            double (*Ma)[4] = M_all[c1];
            double (*Mb)[4] = M_all[c2];
            Nq = generate_next_quaternion(q_sparsity, Q, 1);
            int percent_complete = 5;
            combinations_tried++;
            printf("%d/%d:  0%%  ", combinations_tried, (NM-1)*NM/2);
            fflush(stdout);
            for (int i=0; i<Nq; i++){
                if (100.0*i/Nq >= (double)percent_complete){
                    printf("%d%%  ", percent_complete);
                    fflush(stdout);
                    percent_complete += 5;
                }
                if (generate_next_quaternion(q_sparsity, Q, 0) < 0){
                    printf("There are no more quaternions to search, but Nq still requested more. Exiting.\n");
                    exit(1);
                }
                q = Q;
                inside_S = 1;
                for (int ih=0; ih<NM; ih++){
                    h[ih] = h_func(q, M_all[ih]);
                    if (h[ih] > 0) inside_S = 0;
                }
                if (inside_S && q[3] >= 0){
                    Nw = solve_for_w2(q, Ma, Mb, Nwmin, Nwmax, w_sparsity, W);
                    for (int j=0; j<Nw; j++){
                        w = &W[3*j];
                        Ncon = get_active_constraints(q, w, A, b);
                        if (Ncon == -1) continue; // this means that one of the H values is positive
                        valid_point = is_viable_exact(A,b,Ncon);
                        n_opt++;
                        if (!valid_point){
                            // if (NM==8){
                            //     printf("What is going on? q=[%lf  %lf  %lf  %lf], w = [%lf  %lf  %lf]\n",q[0],q[1],q[2],q[3],w[0],w[1],w[2]);
                            // }
                            memcpy(&track[N_track][0], q, 4*sizeof(double));
                            memcpy(&track[N_track][4], w, 3*sizeof(double));
                            track[N_track][7] = 0.0;
                            for (int ih=0; ih<NM; ih++) track[N_track][7] += h[ih];
                            fprintf(fout, "%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", track[N_track][0], track[N_track][1], track[N_track][2], track[N_track][3], track[N_track][4], track[N_track][5], track[N_track][7], track[N_track][7]);
                            N_track++;

                            if (N_track >= MAX_TRACKED_POINTS){
                                double (*new_points)[8] = malloc(8*(MAX_TRACKED_POINTS+INIT_TRACKED_POINTS)*sizeof(double));
                                memcpy(new_points, track, 8*MAX_TRACKED_POINTS*sizeof(double));
                                MAX_TRACKED_POINTS += INIT_TRACKED_POINTS;
                                free(track);
                                track = new_points;
                                *ptr = track;
                                long total_memory = (MAX_TRACKED_POINTS*8.0 + Nwmax*3)*sizeof(double);
                                printf("More memory required for tracked points. Now using %.4f GB.\n", total_memory/pow(1024, 3.0));
                            }

                            c1 = NM; c2 = NM; // no need to cover all the other cases if this is encountered

                            if (quit > 0){
                                if (quit > 1){
                                    printf("Found an invalid point:\n");
                                    printf("Index Iq = %d\n", i);
                                    printf("q = [%lf; %lf; %lf; %lf];\n",q[0],q[1],q[2],q[3]);
                                    printf("w = [%lf; %lf; %lf];\n",w[0],w[1],w[2]);
                                    for (int iH=0; iH<NM; iH++){
                                        printf("h = %lf,  H = %lf\n", h_func(q,M_all[iH]), H_func(q,w,M_all[iH]));
                                    }
                                    printf("Constraints:\n");
#if 0
                                    for (int iC=0; iC<Ncon; iC++){
                                        for (int k=0; k<3; k++){
                                            if (A[3*iC+k] >=0) printf("+%.8f*u%d ",A[3*iC+k], k+1);
                                            else printf("%.8f*u%d ",A[3*iC+k], k+1);
                                        }
                                        printf("<= %.8f\n",b[iC]);
                                    }
#else
                                    printf("A = [");
                                    for (int iC=0; iC<Ncon; iC++){
                                        for (int k=0; k<3; k++){
                                            printf("%.8f",A[3*iC+k]);
                                            if (k!=2) printf(" ");
                                        }
                                        if (iC != Ncon-1) printf(";\n     ");
                                        else printf("];\n");
                                    }
                                    printf("b = [");
                                    for (int iC=0; iC<Ncon; iC++){
                                        printf("%.8f",b[iC]);
                                        if (iC != Ncon-1) printf("; ");
                                        else printf("];\n");
                                    }
                                    #endif
                                    printf("\n");
                                }
                                return 1;
                            }
                        }
                    }
                    if (Nw > 0) {
                        n_points+=Nw;
                        n_qvalid++;
                    }
                    else n_points++;
                } else {
                    n_points++;
                }
            }       
            printf("100%%\n");
        }
    }

    printf("Generated %d total points, and tested %d boundary points for viability, consisting of %d unique quaternions. Encountered %d unviable points.\n", n_points, n_opt, n_qvalid, N_track);
    fclose(fout);

    free(h);
    return N_track;
}

