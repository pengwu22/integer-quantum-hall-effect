/********************************************************

    FUNCTION! UN-EXECUTABLE ALONE!

    Compute Eigf and Eigv for a single sample of IQHE

 ********************************************************/

// Required Libiaries
// No worries about RE-INCLUDE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // Record Running Time
#include "mkl.h" //LAPACKE_zhpev

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

#ifndef CONSTDEF
#define CONSTDEF 1
    #define PI 3.14159265358979323846
    #define SWITCH_IMPUR 1
#endif

int randsign(void);
int randsign(void){
/*    if((double)rand()/RAND_MAX < 0.5)
        return 1;
    else return -1;*/
    if(genrand64_real2() < 0.5)
        return 1;
    else return -1;
}

void get_impurMat(int impurNum, double L1, double L2, double *impurMatX, double *impurMatY);
void get_impurMat(int impurNum, double L1, double L2, double *impurMatX, double *impurMatY){
    int i;
    for(i = 0; i < impurNum; i++){
        //impurMatX[i] = randuni()*L1;
        impurMatX[i] = genrand64_real2()*L1;
        //impurMatY[i] = randuni()*L2;
        impurMatY[i] = genrand64_real2()*L2;
    }
}

void get_impurInten(int impurNum, double* impurInten);
void get_impurInten(int impurNum, double* impurInten){
    int i;
    if(SWITCH_IMPUR == 1 && impurNum != 0){
        for(i = 0; i < impurNum; i++){
            impurInten[i] = (double) randsign();
        }
    }
    else{
        for(i = 0; i < impurNum; i++){
            impurInten[i] = 0;
        }
    }
}

void math_complexprod(ptr_complex_double a, ptr_complex_double b, ptr_complex_double result);
void math_complexprod(ptr_complex_double a, ptr_complex_double b, ptr_complex_double result){
    result->real = a->real * b->real - a->imag * b->imag;
    result->imag = a->real * b->imag + a->imag * b->real;
}

void math_expimag(double exponent_im, ptr_complex_double result);
void math_expimag(double exponent_im, ptr_complex_double result){
    // The complex result must have had result before this function.
    result->real = cos(exponent_im);
    result->imag = sin(exponent_im);
}

int math_LargeKroDeltaPBC(int r, int s, int Nphi);
int math_LargeKroDeltaPBC(int r, int s, int Nphi){
    int modulus;
    modulus = (r - s) % Nphi;
    if (modulus == 0){
        return 1;
    }
    else return 0;
}

#ifndef PRINTCMP
#define PRINTCMP 1
void print_complex(ptr_complex_double input);
void print_complex(ptr_complex_double input){
    printf("%f\t%f\n", input->real, input->imag);
}

void print_complexvec(int N, ptr_complex_double input);
void print_complexvec(int N, ptr_complex_double input){
    int i;
    for (i = 0; i < N; i++)
        print_complex(&input[i]);
    printf("\n");
}

void print_eigfunction(int N, ptr_complex_double Z, int j);
void print_eigfunction(int N, ptr_complex_double Z, int j){
    int i;
    // Stored in Column
    for(i = 0; i < N; i++){
        print_complex(&Z[(j-1)*N + i]);
    }
    printf("\n");
}
#endif

/** for eigfunc **/
char jobz = 'V'; // Switch 'V': Compute Both Eigv and Eigf
char uplo = 'U'; // Hamiltonian is stored in Upper-tri packed formuplo;


#ifndef THREAD_HAMI_DEF
#define THREAD_HAMI_DEF 1
    //int samplefunction(ptr_complex_double Coeff234, ptr_complex_double Coeff45, int *Coeff6, int Nphi, int quantaConcentration, double impurConcentration, int thetameshcut, int offhead, ptr_complex_double *eigfptrstore, ptr_double *eigvptrstore, double *timing);
    typedef struct samplefunc_params
    {
        int *theta_1;
        int *theta_2;
        int theta_len;

        ptr_complex_double Vmn;
        ptr_complex_double Coeff234;
        ptr_complex_double Coeff45;
        int *Coeff6;

        int Nphi;
        int thetameshcut;
        int offhead;

        ptr_complex_double *eigfptrstore;
        ptr_double *eigvptrstore;

        double *timing;
    } _Hami_params;
    typedef struct samplefunc_params *ptr_Hami_params;
#endif

void *threads_Hami_func(void *params_pre);
// The definition of thread_Hami_func() is placed after samplefuntion().
// But declaration must preceed samplefunction().

int samplefunction(ptr_complex_double Coeff234, ptr_complex_double Coeff45, int *Coeff6, int Nphi, int quantaConcentration, double impurConcentration, int thetameshcut, int offhead, ptr_complex_double *eigfptrstore, ptr_double *eigvptrstore, double *timing, int num_threads, int **thds_theta1, int **thds_theta2, int *theta_len, pthread_t *peer_thds);
int samplefunction(ptr_complex_double Coeff234, ptr_complex_double Coeff45, int *Coeff6, int Nphi, int quantaConcentration, double impurConcentration, int thetameshcut, int offhead, ptr_complex_double *eigfptrstore, ptr_double *eigvptrstore, double *timing, int num_threads, int **thds_theta1, int **thds_theta2, int *theta_len, pthread_t *peer_thds)
{
    // Nphi determines side length
    double Area;
    Area = (double) Nphi / quantaConcentration;

    double L1;
    L1 = sqrt(Area);
    double L2 = L1;

    double thetastep;
    thetastep = 2 * PI / thetameshcut;

    // Matrix: Size N
    int N = Nphi;
    // Loop Index
    int l, k, m, n; // Counting Element
    int i, j; // Reserved Index

/************************* IMPURITY COEFFICIENT **************************/

    int impurNum; // Must integer
    impurNum = (int) floor(impurConcentration * Area + 0.5);

    double *impurMatX;
    double *impurMatY;
    double *impurInten;
    impurMatX = (double *) malloc (impurNum * sizeof(double)); // Must put this in main() instead of any function inside
    impurMatY = (double *) malloc (impurNum * sizeof(double));
    impurInten = (double *) malloc (impurNum * sizeof(double));
    // Random Part
    //srand(seed); INPUT: RANDOM OUTSIDE IN MAIN down to machine
    get_impurMat(impurNum, L1, L2, impurMatX, impurMatY);
    get_impurInten(impurNum, impurInten);


    double CoeffX_im; //imaginary but without i
    double CoeffY_im;
    ptr_complex_double Coeff1buf;
    Coeff1buf = (ptr_complex_double) malloc(sizeof(_complex_double));
    ptr_complex_double buffer_complex;
    buffer_complex = (ptr_complex_double) malloc(sizeof(_complex_double));

    // Prefactor I: Independent of theta1, theta2, l(row), k(column)
    // Potential Vmn (Sum of Complex Exp)

    int dim_n;
    dim_n = offhead+1+offhead;

    ptr_complex_double Vmn;
    Vmn = (ptr_complex_double) malloc ( dim_n * dim_n * sizeof(_complex_double));

    for(m = -offhead; m <= offhead; m++){
        for(n = -offhead; n <= offhead; n++){
        // Prefactor I
            Coeff1buf->real = 0;
            Coeff1buf->imag = 0;
            CoeffX_im = -2*PI*((double)m/L1);// Force double for safty
            CoeffY_im = -2*PI*((double)n/L2);
            for(i = 0; i < impurNum; i++){
                math_expimag(CoeffX_im*impurMatX[i] + CoeffY_im*impurMatY[i], buffer_complex);
                buffer_complex->real *= impurInten[i];// */L1/L2 (Appear in Hami instead)
                buffer_complex->imag *= impurInten[i];// */L1/L2
                Coeff1buf->real += buffer_complex->real;
                Coeff1buf->imag += buffer_complex->imag;
            }
            Vmn[(n+offhead) + dim_n*(m+offhead)] = *Coeff1buf;
        }
    }

/****************************************************************************/


        ptr_Hami_params *paramsptrstore;
        paramsptrstore = (ptr_Hami_params *) malloc ( num_threads * sizeof(ptr_Hami_params) );

        int id = 0;
        for(id = 0; id < num_threads; id++) {
            paramsptrstore[id] = (ptr_Hami_params) malloc(sizeof(_Hami_params));

            paramsptrstore[id]->Vmn = Vmn;
            paramsptrstore[id]->Coeff234 = Coeff234;
            paramsptrstore[id]->Coeff45 = Coeff45;
            paramsptrstore[id]->Coeff6 = Coeff6;

            paramsptrstore[id]->Nphi = Nphi;
            paramsptrstore[id]->thetameshcut = thetameshcut;
            paramsptrstore[id]->offhead = offhead;

            paramsptrstore[id]->eigfptrstore = eigfptrstore;
            paramsptrstore[id]->eigvptrstore = eigvptrstore;

            paramsptrstore[id]->timing = timing;

            paramsptrstore[id]->theta_1 = thds_theta1[id];
            paramsptrstore[id]->theta_2 = thds_theta2[id];
            paramsptrstore[id]->theta_len = theta_len[id];
            pthread_create(&(peer_thds[id]), NULL, threads_Hami_func,  (void* )paramsptrstore[id]);
        }

        // join all the threads
        for(id = 0; id < num_threads; id++)
            pthread_join(peer_thds[id], NULL);

        for(id = 0; id < num_threads; id++)
            free(paramsptrstore[id]);

        free(paramsptrstore);


/****************************************************************************/

    free(impurMatX);
    free(impurMatY);
    free(impurInten);
    free(Coeff1buf);
    free(buffer_complex);
    free(Vmn);

    return 0;
}


void *threads_Hami_func(void *params_pre){

ptr_Hami_params params = (ptr_Hami_params) params_pre;

    // Program running time
    clock_t tic;
    clock_t toc;
        // Main Parameters: INPUT
    //int thetameshcut = 10; // !! Var4
    int thetameshcut;
    thetameshcut = params->thetameshcut;
    double thetastep;
    thetastep = 2 * PI / thetameshcut;
    //int offhead = 10; // Var5
    int offhead;
    offhead = params->offhead;

    // Matrix Part: Size N
    int N = params->Nphi;

    // Loop Index
    int l, k, m, n; // Counting Element
    int i, j; // Reserved Index


    int dim_n, dim_nm, dim_ml, dim_nm1, dim_ml2;
    dim_n = offhead+1+offhead;
    dim_nm = dim_n * dim_n;
    dim_nm1 = dim_nm * (thetameshcut + 1);
    dim_ml = dim_n * N;
    dim_ml2 = dim_ml * (thetameshcut + 1);

    ptr_complex_double Hamiltonian; // Hamiltonian AP // Defined as input argument
    Hamiltonian = (ptr_complex_double) malloc (N * (N+1) / 2 * sizeof(_complex_double));

    ptr_complex_double mnSum;
    mnSum = (ptr_complex_double) malloc(sizeof(_complex_double));

    ptr_complex_double buffer_complex_pre;
    buffer_complex_pre = (ptr_complex_double) malloc(sizeof(_complex_double));
    ptr_complex_double buffer_complex;
    buffer_complex = (ptr_complex_double) malloc(sizeof(_complex_double));
    int Coeff6buf;
    int index_Hami = 0;
    int index_nm = 0;


/************************/
    int counttheta1 = 0;
    int counttheta2 = 0;

    int lenthread = 0;
    for(lenthread = 0; lenthread < params->theta_len; lenthread++) {
        counttheta1 = params->theta_1[lenthread];
        counttheta2 = params->theta_2[lenthread];

        /********************************/
        /* The Most Time-Consuming Part */
        /********************************/

        tic = clock();// Timing Start: Sec1
        for(l = 1; l <= N; l++){
            for(k = l; k <= N; k++){
                index_Hami = l+(k-1)*k/2-1; //Upper Triangular Storage
                mnSum->real = 0;
                mnSum->imag = 0;

                index_nm = 0;
                for(m = -offhead; m <= offhead; m++){
                    for(n = -offhead; n <= offhead; n++){

                        Coeff6buf = params->Coeff6[n+offhead+dim_n*index_Hami];
                        if(Coeff6buf == 0){
                            index_nm++;
                            continue;
                        }

                        //Coeff234[(n+offhead) + dim_n*(m+offhead) + dim_nm*counttheta1] = *buffer_complex;
                        //Coeff45[(m+offhead) + dim_n*(l-1) + dim_ml*counttheta2] = *buffer_complex;
                        math_complexprod(&(params->Coeff234[index_nm + dim_nm*counttheta1]), &(params->Vmn[index_nm]), buffer_complex_pre);
                        math_complexprod(buffer_complex_pre, &(params->Coeff45[(m+offhead) + dim_n*(l-1) + dim_ml*counttheta2]), buffer_complex);
                        index_nm++;

                        // C1*C2*C3*C4*C5*C6 = C12345*C6
                        mnSum->real += buffer_complex->real * Coeff6buf;
                        mnSum->imag += buffer_complex->imag * Coeff6buf;
                    }
                }

                Hamiltonian[index_Hami] = *mnSum;
            }
        }
        toc = clock(); // Timing end: Sec1
        params->timing[1] += (double)(toc - tic) / CLOCKS_PER_SEC;

        // Store eigf and eigv
        params->eigfptrstore[counttheta1 * (thetameshcut + 1) + counttheta2] = (ptr_complex_double) malloc ( N * N * sizeof(_complex_double));
        params->eigvptrstore[counttheta1 * (thetameshcut + 1) + counttheta2] = (double *) malloc (N * sizeof(double));

        // Diagonalization
        tic = clock(); // Timing Start: Sec2
        LAPACKE_zhpev(LAPACK_COL_MAJOR, jobz, uplo, N, Hamiltonian, params->eigvptrstore[counttheta1 * (thetameshcut + 1) + counttheta2], params->eigfptrstore[counttheta1 * (thetameshcut + 1) + counttheta2], N);
        toc = clock(); // Timing End: Sec2
        params->timing[2] += (double)(toc - tic) / CLOCKS_PER_SEC;
    }


    free(Hamiltonian);
    free(mnSum);
    free(buffer_complex_pre);
    free(buffer_complex);
    // DO NOT f_r_e_e the pointers stored in eigf/vstore.

    //printf("I am a nice thread.\n");

    pthread_exit((void *) 0);
    return (void*) 0;// For compiler syntax check
}
