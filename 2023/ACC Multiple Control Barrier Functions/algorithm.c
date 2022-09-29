#include "verify.h"

void construct_new_M(double v[3], double M[4][4], double theta){
    memset(M, 0, 16*sizeof(double));

    double magr = Mag3(rb);
    double magv = Mag3(rb);
    if (fabs(magr) < eps_tol || fabs(magv) < eps_tol){
        printf("Vectors in M matrix are zero magnitude. Exiting.\n");
        exit(1);
    }
    for (int i=0; i<3; i++){
        rb[i] /= magr;
        v[i] /= magv;
    }

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            M[i][j] += v[i]*rb[j] + rb[i]*v[j];
        }
        M[i][i] -= Dot3(v, rb) + cos(theta);
    }

    double temp[3];
    Cross3(rb, v, temp);
    for (int i=0; i<3; i++){
        M[3][i] = temp[i];
        M[i][3] = temp[i];
    }
    M[3][3] = Dot3(rb, v) - cos(theta);
}

int get_cluster(double (*input)[8], int num_points, double (*output)[8]){
    int count = 0;
    char filename[20];
    sprintf(filename, "data/cluster%d.dat", NM);
    FILE *fout = fopen(filename, "w+");

    memcpy(output[count], input[0], 8*sizeof(double));
    fprintf(fout, "%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", input[0][0], input[0][1], input[0][2], input[0][3], input[0][4], input[0][5], input[0][7], input[0][7]);
    memset(input[0], 0, 8*sizeof(double));
    count++;
    for (int k=0; k<K_MAX_CLUSTER; k++){
        int points_to_compare = count;
        for (int i=1; i<num_points; i++){
            double v1[3];
            RotQ(rb, input[i], v1);
            for (int j=0; j<points_to_compare; j++){
                if (input[i][3] != 0){
#ifdef _RUN_CASE_3_
                    /**
                     * The following lines result in covering about half of the total quaternions 
                     * because there exists an extra degree of freedom in the quaternions compared to the pointing vectors.
                     * The algorithm will work just as well with either code, but the latter section is more intuitive.
                     * */
                    double q_diff[4];
                    QxQT(input[i], output[j], q_diff);
                    if (1.0 - fabs(q_diff[3]) < q_neighbor_tol){
#else
                    double v2[3], v_diff[3];
                    RotQ(rb, output[j], v2);
                    for (int p=0; p<3; p++) v_diff[p] = v1[p] - v2[p];
                    if (Mag3(v_diff) < q_neighbor_tol){
#endif
                        memcpy(output[count], input[i], 8*sizeof(double));
                        fprintf(fout, "%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", input[i][0], input[i][1], input[i][2], input[i][3], input[i][4], input[i][5], input[i][6], input[i][7]);
                        memset(input[i], 0, 8*sizeof(double));
                        count++;
                        break;
                    }
                }
            }
        }
        if (count==points_to_compare) {
            k = K_MAX_CLUSTER;
            printf("Cluster successfully found!\n");
        }
    }
    printf("Returned cluster of size %d\n", count);
    fclose(fout);
    return count;
}

void get_cbf(double (*cluster)[8], int num_points, double out[4][4]){
    static FILE *fout = NULL;
    int min_index = 0;
    int max_index = 0;

    for (int i=0; i<num_points; i++){
        if (cluster[i][7] > cluster[min_index][7]){ // compare the h values and look for the point most on the boundary
            min_index = i;
        }
        if (cluster[i][7] < cluster[min_index][7]){ // compare the h values and look for the point most interior
            max_index = i;
        }
    }

    double q_bd[4], q_int[4];
    for (int i=0; i<num_points; i++){
        for (int j=0; j<4; j++){
            q_int[j] += ((cluster[i][7] - cluster[min_index][7])/cluster[max_index][7])*cluster[i][j];
            q_bd[j] += (1.0 - (cluster[i][7] - cluster[min_index][7])/cluster[max_index][7])*cluster[i][j];
        }
    }
    double mag = sqrt(q_int[0]*q_int[0] + q_int[1]*q_int[1] + q_int[2]*q_int[2] + q_int[3]*q_int[3]);
    for (int i=0; i<4; i++) q_int[i] /= mag;
    mag = sqrt(q_bd[0]*q_bd[0] + q_bd[1]*q_bd[1] + q_bd[2]*q_bd[2] + q_bd[3]*q_bd[3]);
    for (int i=0; i<4; i++) q_bd[i] /= mag;

    int num_not_in_cluster = 1;
    double theta = 2*acos(-1)/180; // guess that 2 degrees is a good solution
    double v_remove[3], v_bd[3], v_int[3];
    double M_new[4][4];
#ifdef _RUN_CASE_1_
    /* Either these two anchor points, or the anchor points used in _RUN_CASE_2_ will work in this case. 
     * We choose this particular case just to demonstrate that there are system-specific options when designing getCBF
     */
    RotQ(rb, cluster[min_index], v_bd);
    RotQ(rb, cluster[max_index], v_int);
#endif
#ifdef _RUN_CASE_2_
    RotQ(rb, q_bd, v_bd);
    RotQ(rb, q_int, v_int);
#endif
#ifdef _RUN_CASE_3_
    /* Because there is only one cluster in this case, the "most interior" point could be on either side, so we use a fixed anchor point instead */
    RotQ(rb, cluster[min_index], v_bd);
    v_int[0] = 0; v_int[1] = 0; v_int[2] = v_bd[2]/fabs(v_bd[2]);
#endif
    printf("Initial vectors are v_remove = [%lf  %lf  %lf] and v_opp = [%lf  %lf  %lf]\n", v_bd[0], v_bd[1], v_bd[2], v_int[0], v_int[1], v_int[2]);
    memcpy(v_remove, v_bd, 3*sizeof(double));
    while (num_not_in_cluster > 0){
        construct_new_M(v_remove, M_new, theta);
        num_not_in_cluster = 0;
        for (int i=0; i<num_points; i++){
            if (h_func(&cluster[i][0],M_new) <= 0 && H_func(&cluster[i][0],&cluster[i][4],M_new) <= 0){
                num_not_in_cluster++;
                break;
            }
        }
        
        if (num_not_in_cluster > 0){
            double delta = 1*acos(-1)/180;
            double axis[3], mag, delta_q[4], new_v[3];
            theta += delta;
            Cross3(v_int, v_bd, axis);
            mag = Mag3(axis);
            if (mag < eps_tol) printf("Axis of rotation is zero in get_cbf.\n");
            for (int k=0; k<3; k++) delta_q[k] = axis[k]/mag*sin(delta/2);
            delta_q[3] = cos(delta/2);
            RotQ(v_remove, delta_q, new_v);
            memcpy(v_remove, new_v, 3*sizeof(double));
            printf("Performing another iteration in get_cbf. Current theta = %lf\n", theta);
        }
    }

#ifdef _RUN_CASE_2_
    /* This section checks if the new CBF will add a conflict with any of the already existing CBFs.
     * Without this extra step, _RUN_CASE_2_ will take considerably longer.
     * This section of code is essentially an instance of hyperparameter tuning. */
    num_not_in_cluster = 1;
    while (num_not_in_cluster > 0){
        construct_new_M(v_remove, M_new, theta);

        int NM_temp = NM;
        double (*M_temp)[4][4] = malloc(16*NM_temp*sizeof(double));
        memcpy(M_temp, M_all, 16*NM_temp*sizeof(double));
        double (*tracked)[8] = malloc(8*MAX_TRACKED_POINTS*sizeof(double));

        printf("Checking %d combinations for possible future conflicts.\n", NM_temp);
        for (int c=0; c<NM_temp; c++){
            memcpy(M_all[0], M_temp[c], 16*sizeof(double));
            memcpy(M_all[1], M_new, 16*sizeof(double));
            NM = 2;
            num_not_in_cluster = check_points(&tracked, 1);
            if (num_not_in_cluster > 0) break;
        }
        free(tracked);
        memcpy(M_all, M_temp, 16*NM_temp*sizeof(double));
        NM = NM_temp;

        if (num_not_in_cluster > 0){
            double delta = 1*acos(-1)/180;
            double axis[3], mag, delta_q[4], new_v[3];
            theta += delta;
            Cross3(v_int, v_bd, axis);
            mag = Mag3(axis);
            if (mag < eps_tol) printf("Axis of rotation is zero in get_cbf.\n");
            for (int k=0; k<3; k++) delta_q[k] = axis[k]/mag*sin(delta/2);
            delta_q[3] = cos(delta/2);
            RotQ(v_remove, delta_q, new_v);
            memcpy(v_remove, new_v, 3*sizeof(double));
            printf("Performing another iteration in get_cbf. Current theta = %lf\n", theta);
        }

    }
#endif
    

    memcpy(out, M_new, 16*sizeof(double));
    
    printf("Added constraint #%d at v = [%lf  %lf  %lf] and theta = %lf\n", NM+1, v_remove[0], v_remove[1], v_remove[2], theta);
    if (fout==NULL) fout = fopen("data/vectors.dat", "w+");
    fprintf(fout, "%lf %lf %lf %lf\n", v_remove[0], v_remove[1], v_remove[2], theta);
    fflush(fout);
}

void print_solution(void){
    FILE *fout = fopen("data/solution.dat", "w+");
    fprintf(fout, "double M_all[%d][4][4] = {{{", NM);
    for (int m=0; m<NM; m++){
        for (int i=0; i<4; i++){
            for (int j=0; j<4; j++){
                fprintf(fout, "%.16f", M_all[m][i][j]);
                if (j<3) fprintf(fout, ", ");
            }
            if (i<3) fprintf(fout, "},\n                                    {");
        }
        if (m<NM-1) fprintf(fout, "}},\n                                   {{");
        else fprintf(fout, "}}};\n");
    }
    fclose(fout);
}

int main(){
    clock_t start, end;
     
    start = clock();

    MAX_TRACKED_POINTS = INIT_TRACKED_POINTS;
    double (*tracked)[8] = malloc(8*MAX_TRACKED_POINTS*sizeof(double));
    M_all = malloc(16*NM0*sizeof(double));
    memcpy(M_all, M_start, 16*NM0*sizeof(double));
    NM = NM0;

    int N_conflict;

    char all_good = 0;
    while (!all_good){
        N_conflict = check_points(&tracked, 0);
        
        if (N_conflict==0){
            all_good = 1;
            break;
        }

        // Find where we are going to put the next constraint
        double (*M_new)[4][4] = malloc((NM+1)*16*sizeof(double));
        memcpy(M_new, M_all, NM*16*sizeof(double));
        free(M_all);
        M_all = M_new;

        double (*cluster)[8] = malloc(8*N_conflict*sizeof(double));
        int N_cluster = get_cluster(tracked, N_conflict, cluster);
        get_cbf(cluster, N_cluster, M_new[NM]);
        free(cluster);
        NM += 1;

    }

    print_solution();


    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total Program Run Time was %lf seconds.\n", cpu_time_used);
    return 1;
}