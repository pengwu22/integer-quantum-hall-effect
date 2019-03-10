/********************************************************

    MAIN PROGRAM:

    TITLE: Integer Hall Pleateau Transition
    AUTHOR: Peng Wu @ Zhejiang University Physics Dprt.
    DATED: Dec. 29, 2015

    To compile:
    gcc main.c -o main -mkl=sequential -lpthread

    P.S.
    -mkl=sequential: Intel MKL sequential computing
    -lpthread: Multi-thread package
    -llapack: Lapack package
    -lblas: Blas package
    -lm: math.h

 ********************************************************/

// Every <*.h> here is Stardard C Library Header.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "mt64.h"
#include "mkl.h" // -mkl
#include "pthread.h" // -lpthread

// Every "*.c" here is a section of script with static order.
#include "mt64.c" // included before samplefunc.c
#include "samplefunc_threads.c" // included before presamplefunc.c
#include "presamplefunc.c"
#include "chernfunc.c"
#include "argv_enable.c"

/***** READ SEED OR NOT? ****/
#define SWTICH_READSEED 0
// 0: Auto generate random seed // DEFAULT
// 1: Read seed from file "seed"
/****************************/

/* INPUT ARGVS FROM OUTSIDE? */
#define SWITCH_MAINARGV 1 // Switch of main() argument input
// 1: Must input from outside // DEFAULT
// 0: No argument input (for "gdb")
/****************************/

/******* TOTAL THREADS ******/
#define NUM_THREADS 4
/****************************/

/************************************************************/
/****** Think twice before modify the following part! *******/
/************************************************************/

#ifndef STRUCTDEF
#define STRUCTDEF 1
    typedef MKL_Complex16 _complex_double;
    typedef MKL_Complex16 *ptr_complex_double;
    //typedef struct cmplxnode *ptr_complex_double;
    //typedef struct cmplxnode { double real, imag; } _complex_double;
#endif

#ifndef PTRDBDEF
#define PTRDBDEF 1
    typedef double *ptr_double;
#endif

#define MAX_SAMPLENUM 9999999
#define MAX_TIMEINGSEC 4

int main(int argc, char *argv[])
{
    FILE *output;
    FILE *fp_seed;

/************************************************************/
/********************* INITIALIZATION ***********************/
/************************************************************/

    // Main Variables
    // Var1: Number of Flux Quanta
    int Nphi = 8;
    // Var2: Concentration of Flux Quanta
    int quantaConcentration = 1; // Fixed = 1
    // Var3: Concentration of Impurities
    double impurConcentration = 16; // C_impur > C_quanta
    // Var4: Side Partition of Meshgrid in Theta Space
    int thetameshcut = 10;
    // Var5: Fourier Series Summation Cutoff Head
        // Old Version: // int offhead = 7;
        // New Version:
    int MoverL = 5; // ~o( sqrt(Cimp) ), but never less than sqrt(Cimp), Larger the Better
    int offhead = 10; // init. with a reasonable value to protect the program
    offhead = ceil( MoverL * sqrt(Nphi) / 2 ); // Revalued
    // Var6: Number of Samples
    int sampletotal = 1;

    if (SWITCH_MAINARGV == 1){
        int impurConcentrationUp;
        int impurConcentrationDw;
        Nphi = mystr2int(argv[1]); // Var1: Must Integer
        quantaConcentration = 1; // Var2: Must Integer
        impurConcentrationUp = mystr2int(argv[2]); // Var3: Double C_impur > C_quanta
        impurConcentrationDw = mystr2int(argv[3]); // Var3: Double C_impur > C_quanta
            impurConcentration = (double)impurConcentrationUp; // impurConcentrationDw;
        thetameshcut = mystr2int(argv[4]); // !! Var4
        //offhead = mystr2int(argv[5]); // Var5
        //MoverL = mystr2int(argv[5]); // Var5
        MoverL = ceil( sqrt(impurConcentration) ) + 1; // Var5: Fourier Series Precision
            offhead = ceil( MoverL * sqrt(Nphi) / 2 ); // Revalued
        sampletotal = mystr2int(argv[5]);
    }

    // Variable Initialization Error check
    if(sampletotal > MAX_SAMPLENUM){
        printf("Too Many Sample!\n");
        printf("The number of sample should be less that %d\n", MAX_SAMPLENUM);
        return 0;
    }

    // Show on screen
    printf("\nPROGRAM IS INITIALIZING...\n");
    printf("<int> Nphi = %d\n", Nphi);
    //printf("<int> Quanta Concentration = %d\n", quantaConcentration);
    printf("<double> Impur Concentration = %f\n", impurConcentration);
    printf("<int> Theta-Space One-Side Meshgrid = %d\n", thetameshcut);
    //printf("<int> M/L = %d\t(Must Exceed %f)\n", MoverL, sqrt(impurConcentration));
    printf("<int> Sample# = %d\n", sampletotal);
    printf("NOTE: If any argument above is not the one you input, please stop the program immediately!\n\n");

    printf("CALCULATION STARTS NOW!\n");

/************************************************************/

    // MAIN STORAGE
    ptr_complex_double *eigfptrstore;
    eigfptrstore = (ptr_complex_double *) malloc ((thetameshcut+1) * (thetameshcut+1) * sizeof(ptr_complex_double));
    ptr_double *eigvptrstore;
    eigvptrstore = (ptr_double *) malloc ((thetameshcut+1) * (thetameshcut+1) * sizeof(ptr_double));

    // Record Program Running Time
    clock_t tic, toc;
    double *timing;
    timing = (double *) malloc (MAX_TIMEINGSEC * sizeof(double));
    timing[0] = 0;
    timing[1] = 0;
    timing[2] = 0;
    timing[3] = 0;

    // Initialize a random seed
    unsigned int seed;
    unsigned int offset = 0; // NEW
    if (SWTICH_READSEED == 1){
        int status_fscanf;
        fp_seed = fopen("./seed","r");
        status_fscanf = fscanf(fp_seed,"%u",&seed);
        if(status_fscanf == EOF){
            printf("Seed Scaning Error.\n");
            return 0;
        }
        printf("Read-in SEED = %u\n\n", seed);
        fclose(fp_seed);
    }
    else{
        srand((unsigned)time(NULL));
        seed = rand(); //ATTENTION FOR GROUP SUBMITTING
        printf("Pseudo-random seed is generated upon time.\n");
    }

    // mt19937-64 seed initialization
    init_genrand64(seed*113ULL+13127641ULL*offset);

    // Save the initialized random seed
    output = fopen("seed", "w");
    fprintf(output, "%u\n", seed);
    fclose(output);

    // Ready to record error-chern sample
    int *howmanyerrchern;
    howmanyerrchern = (int *) malloc (sizeof(int));
    *howmanyerrchern = 0;

    // Ready to count chern
    int *howmanyext;
    howmanyext = (int *) malloc (sizeof(int));
    *howmanyext = 0;

    // Renew New Output File for A Group of Samples
    output = fopen ("ensemble_data","w");
    fclose(output);

/************************************************************/
/********************* PREPARATION  *************************/
/************************************************************/
    int dim_n, dim_nm, dim_ml, dim_nm1, dim_ml2;
    dim_n = offhead+1+offhead;
    dim_nm = dim_n * dim_n;
    dim_nm1 = dim_nm * (thetameshcut + 1);
    dim_ml = dim_n * Nphi;
    dim_ml2 = dim_ml * (thetameshcut + 1);

    ptr_complex_double Coeff234; //_nmtheta1
    ptr_complex_double Coeff45; //_mltheta2
    Coeff234 = (ptr_complex_double) malloc (dim_nm1 * sizeof(_complex_double));
    Coeff45 = (ptr_complex_double) malloc (dim_ml2 * sizeof(_complex_double));

    int *Coeff6;
    Coeff6 = (int *) malloc (dim_n * Nphi * (Nphi+1) / 2 * sizeof(int));
/************************************************************/
    printf("Preparing...\n");
    tic = clock();

    presamplefunc(Coeff234, Coeff45, Coeff6, Nphi, quantaConcentration, impurConcentration, thetameshcut, offhead);

    toc = clock();
    timing[0] = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("Preparation Ends.\n");
/************************************************************/
/********************* MAIN PART ****************************/
/************************************************************/

    int i, j;

/******* Split theta space to thread *******/

    int num_threads = NUM_THREADS;
    // paramsters for peer_solve
    int **thds_theta1, **thds_theta2, *theta_len;
    thds_theta1 = (int**) malloc (num_threads * sizeof(int *));
    thds_theta2 = (int**) malloc (num_threads * sizeof(int *));

    for(i = 0; i < num_threads; i++) {
        thds_theta1[i] = (int *) malloc ( ( (thetameshcut + 1) * (thetameshcut + 1) / num_threads + 1 ) * sizeof(int) );
        thds_theta2[i] = (int *) malloc ( ( (thetameshcut + 1) * (thetameshcut + 1) / num_threads + 1 ) * sizeof(int) );
    }
    theta_len = (int*) malloc (num_threads * sizeof(int));

    memset(theta_len, 0, num_threads*sizeof(int));
    for(i = 0; i < (thetameshcut + 1); i++)
        for(j = 0; j < (thetameshcut + 1); j++) {
            thds_theta1[(i * (thetameshcut + 1) + j) % num_threads][theta_len[(i * (thetameshcut + 1) + j) % num_threads]] = i;
            thds_theta2[(i * (thetameshcut + 1) + j) % num_threads][theta_len[(i * (thetameshcut + 1) + j) % num_threads]++] = j;
        }

/********* !! MAIN LOOP !! ******************/

    pthread_t *peer_thds;
    peer_thds =  (pthread_t *) malloc (num_threads*sizeof(pthread_t));

    for(i = 1; i <= sampletotal; i++){
        // Show Where I am
        printf("SAMPLE%8d\n", i);

        // Generate Eigf
        samplefunction(Coeff234, Coeff45, Coeff6, Nphi, quantaConcentration, impurConcentration, thetameshcut, offhead, eigfptrstore, eigvptrstore, timing, num_threads, thds_theta1, thds_theta2, theta_len, peer_thds);

        // Calculate Chern
        tic = clock();
        chernfunc(i, Nphi, thetameshcut, eigfptrstore, eigvptrstore, howmanyerrchern, howmanyext);
        toc = clock();
        timing[3] += (double)(toc - tic) / CLOCKS_PER_SEC;
    }

    free(peer_thds);

/************************************************************/

    for(i = 0; i < num_threads; i++) {
        free(thds_theta1[i]);
        free(thds_theta2[i]);
    }
    free(thds_theta1);
    free(thds_theta2);
    free(theta_len);

/************************************************************/

    // Write ensemble_info
    output = fopen("ensemble_info","w");
    fprintf(output, "%6d%9.1f%6d%6d%9d%6d%10d%12.2f%12.2f%12.2f%12.2f\n", Nphi, impurConcentration, thetameshcut, MoverL, sampletotal, *howmanyerrchern, *howmanyext, timing[0], timing[1], timing[2], timing[3]);
    fclose(output);

    // Update on screen
    printf("\n");
    printf("RESULTS:\n");
    //printf("Recheck arguments, and if any changed, abandon the results: \n");
    //printf("%6d%8.3f%6d%6d%8d\n\n", Nphi, impurConcentration, thetameshcut, offhead, sampletotal);
    printf("Sample# with error total conductance / Total Sample# = %d / %d.\n", *howmanyerrchern, sampletotal);
    printf("Extended States / All states = %d / %d.\n", *howmanyext, sampletotal*Nphi);
    printf("Prepa_time =%12.2f\nHamit_time =%12.2f\nDiagn_time =%12.2f\nChern_time =%12.2f\n", timing[0], timing[1]/num_threads, timing[2]/num_threads, timing[3]);
    printf("\nPROGRAM ENDS.\n");

    // Empty Space
    free(eigfptrstore);
    free(eigvptrstore);
    free(timing);
    free(howmanyerrchern);
    free(howmanyext);
    free(Coeff234);
    free(Coeff45);
    free(Coeff6);

    return 0;
}
