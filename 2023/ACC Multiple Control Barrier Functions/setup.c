#include "verify.h"

#define M_INDEX (0)

/* This program checks whether a single constraint is a CBF. If not, then we have to adjust the value of mu or w_max */
int main(){
    clock_t start, end;
     
    start = clock();

    double (*M)[4] = M_start[M_INDEX];
    double *Q;
    q_sparsity = 0.1;
    int Nq = generate_unit_quaternions(q_sparsity, &Q);

    int Nmax = 10000, Nw;
    w_sparsity = 2.0*w_max/100;
    double *W = malloc(3*Nmax*sizeof(double));
    double h, *q, *w, u[3], un, Hdot;
    int count = 0;
    int count_w = 0;
    int count_q = 0;

    int percent_complete = 5;
    for (int i=0; i<Nq; i++){
        if (100.0*i/Nq >= (double)percent_complete){
            printf("%d%%  ", percent_complete);
            fflush(stdout);
            percent_complete += 5;
        }
        q = &Q[4*i];
        h = h_func(q, M);
        if (h <= 0){
            count_q++;
            Nw = generate_w1(q, M, w_sparsity, Nmax, W);
            for (int j=0; j<Nw; j++){
                count_w++;
                w = &W[3*j];
                hdot_affine(q, M, u);
                un = Mag3(u);
                if (un > eps_tol){
                    for (int p=0; p<3; p++) u[p] *= -u_max/(un*sqrt(3)); 
                    Hdot = Hdot_func(q, w, u, M);
                    if (Hdot > 0){
                        printf("Found an invalid point:\n");
                        printf("q = [%lf; %lf; %lf; %lf];\n",q[0],q[1],q[2],q[3]);
                        printf("w = [%lf; %lf; %lf];\n",w[0],w[1],w[2]);
                        printf("h = %lf,  H = %lf\n", h, H_func(q,w,M));
                        printf("Hdot = %lf\n", Hdot);
                        printf("\n");
                        exit(1);
                    }
                } else {
                    printf("u_guess = 0. You should look into that. Exiting.\n");
                    exit(1);
                }
            }
            if (Nw > 0) count+=Nw;
            else count++;
        } else {
            count++;
        }
    }
    printf("%d%%\n\n", percent_complete);




    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total Program Run Time was %lf seconds.\n", cpu_time_used);
    printf("Generated %d quaternions, %d total points, and tested %d boundary points for viability.\n", Nq, count, count_w);
    printf("Tested %d unique quaternion points.\n", count_q);
    return 1;
}

