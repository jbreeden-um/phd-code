#include "verify.h"

/* This is a by-hand solution based on visualizing the problem */
#define N_expert (1)
double M_expert[N_expert][4][4] = {{{-0.898794046299167,   1.000000000000000,                   0,                   0},
                                    { 1.000000000000000,  -0.898794046299167,                   0,                   0},
                                    {                 0,                   0,  -1.898794046299167,                   0},
                                    {                 0,                   0,                   0,   0.101205953700833}}};

/* This is the output of _RUN_CASE_2_ */
// #define N_expert (2)
// double M_expert[N_expert][4][4] = {{{-0.9038894301175082, 0.9985825983165566, -0.0375961825799410, -0.0375961825799410},
//                                     {0.9985825983165566, -0.9087261439557915, -0.0375961825799410, 0.0375961825799410},
//                                     {-0.0375961825799410, -0.0375961825799410, -1.9048903853532064, -0.0024183569191417},
//                                     {-0.0375961825799410, 0.0375961825799410, -0.0024183569191417, 0.0922748112799068}},
//                                    {{-0.8982391248763710, 0.9984315976012965, 0.0391742364460889, 0.0391742364460889},
//                                     {0.9984315976012965, -0.9143764491969286, 0.0391742364460889, -0.0391742364460889},
//                                     {0.0391742364460889, 0.0391742364460889, -1.9047393846379463, -0.0080686621602788},
//                                     {0.0391742364460889, -0.0391742364460889, -0.0080686621602788, 0.0921238105646467}}};

/* This is the output of _RUN_CASE_3_ */
// #define N_expert (1)
// double M_expert[N_expert][4][4] = {{{-0.8915065534318422, 0.9999870541991667, 0.0035806008456685, 0.0035806008456685},
//                                     {0.9999870541991667, -0.8905064949448935, 0.0035806008456685, -0.0035806008456685},
//                                     {0.0035806008456685, 0.0035806008456685, -1.8909935783875345, 0.0005000292434743},
//                                     {0.0035806008456685, -0.0035806008456685, 0.0005000292434743, 0.1089805300107989}}};

/* This program checks whether a set of pre-generated CBFs indeed form a viability domain */
int main(){
    clock_t start, end;
     
    start = clock();

    MAX_TRACKED_POINTS = INIT_TRACKED_POINTS;
    double (*tracked)[8] = malloc(8*MAX_TRACKED_POINTS*sizeof(double));
    
    NM = 2+N_expert;
    M_all = malloc(16*NM*sizeof(double));
    memcpy(M_all, M_start, 16*2*sizeof(double));
    memcpy(M_all[2], M_expert, 16*N_expert*sizeof(double));

    int N_conflict;

    N_conflict = check_points(&tracked, 2);

    if (N_conflict==0) printf("Successfully checked the proposed solution, and found no conflicts.\n");
    else printf("Conflicts were found. The proposed solution does not yield a viability domain.\n");

    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total Program Run Time was %lf seconds.\n", cpu_time_used);
    return 1;
}