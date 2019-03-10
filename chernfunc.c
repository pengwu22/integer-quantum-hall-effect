#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

void math_conjdotprod(int N, ptr_complex_double x, ptr_complex_double y, ptr_complex_double result);
void math_conjdotprod(int N, ptr_complex_double x, ptr_complex_double y, ptr_complex_double result){
    int i;
    // Refresh first
    result->real = 0;
    result->imag = 0;
    for(i = 0; i < N; i++){
        result->real += x[i].real * y[i].real + x[i].imag * y[i].imag;
        result->imag += x[i].real * y[i].imag - x[i].imag * y[i].real;
    }
}

void readeigf(int N, ptr_complex_double buffer_umat, int m, ptr_complex_double buffer_eigf);
void readeigf(int N, ptr_complex_double buffer_umat, int m, ptr_complex_double buffer_eigf){
    int i;
    for(i = 0; i < N; i++){
        buffer_eigf[i] = buffer_umat[m * N + i];
    }
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

// MAIN FUNCTION
int chernfunc(int samplenum, int Nphi, int thetameshcut, ptr_complex_double *eigfptrstore, ptr_double *eigvptrstore, int *howmanyerrchern, int *howmanyext);
int chernfunc(int samplenum, int Nphi, int thetameshcut, ptr_complex_double *eigfptrstore, ptr_double *eigvptrstore, int *howmanyerrchern, int *howmanyext)
{
    int i,j;

    // Main Var
    double thetastep; // Var4.1 determined by 4
    thetastep = 2 * PI / thetameshcut;

    // Matrix Part: Size N
    int N = Nphi;

/************************************************************/
/***************  MAIN PART  ********************************/
/************************************************************/
    // Initialize Chern Storage
    int m; //State#
    double gammaelement;
    double *gamma;
    gamma = (double *) malloc (Nphi * sizeof(double));

    // Initialize Loop Math Storage
    ptr_complex_double buffer_umat;

    ptr_complex_double buffer_eigf00;
    buffer_eigf00 = (ptr_complex_double) malloc (N * sizeof(_complex_double));
    ptr_complex_double buffer_eigf10;
    buffer_eigf10 = (ptr_complex_double) malloc (N * sizeof(_complex_double));
    ptr_complex_double buffer_eigf01;
    buffer_eigf01 = (ptr_complex_double) malloc (N * sizeof(_complex_double));
    ptr_complex_double buffer_eigf11;
    buffer_eigf11 = (ptr_complex_double) malloc (N * sizeof(_complex_double));

    ptr_complex_double lap1;
    lap1 = (ptr_complex_double) malloc (sizeof(_complex_double));
    ptr_complex_double lap2;
    lap2 = (ptr_complex_double) malloc (sizeof(_complex_double));
    ptr_complex_double lap3;
    lap3 = (ptr_complex_double) malloc (sizeof(_complex_double));
    ptr_complex_double lap4;
    lap4 = (ptr_complex_double) malloc (sizeof(_complex_double));

    double ImLoglap1, ImLoglap2, ImLoglap3, ImLoglap4;

    // loop
    int counttheta1 = 0;
    int counttheta2 = 0;
    for(m = 0; m < N; m++){
        gamma[m] = 0;
        for(counttheta1 = 0; counttheta1 < thetameshcut; counttheta1++){
            for(counttheta2 = 0; counttheta2 < thetameshcut; counttheta2++){

                //Load eigf(0,0) to buffer
                buffer_umat = eigfptrstore[counttheta1 * (thetameshcut + 1) + counttheta2];
                readeigf(N, buffer_umat, m, buffer_eigf00);
                //Load eigf(1,0) to buffer
                buffer_umat = eigfptrstore[(counttheta1 + 1) * (thetameshcut + 1) + counttheta2];
                readeigf(N, buffer_umat, m, buffer_eigf10);
                //Load eigf(0,1) to buffer
                buffer_umat = eigfptrstore[counttheta1 * (thetameshcut + 1) + counttheta2 + 1];
                readeigf(N, buffer_umat, m, buffer_eigf01);
                //Load eigf(1,1) to buffer
                buffer_umat = eigfptrstore[(counttheta1 + 1) * (thetameshcut + 1) + counttheta2 + 1];
                readeigf(N, buffer_umat, m, buffer_eigf11);

                math_conjdotprod(N, buffer_eigf00, buffer_eigf10, lap1);//inside allocate room
                math_conjdotprod(N, buffer_eigf10, buffer_eigf11, lap2);
                math_conjdotprod(N, buffer_eigf11, buffer_eigf01, lap3);
                math_conjdotprod(N, buffer_eigf01, buffer_eigf00, lap4);

                ImLoglap1 = atan2(lap1->imag, lap1->real);
                ImLoglap2 = atan2(lap2->imag, lap2->real);
                ImLoglap3 = atan2(lap3->imag, lap3->real);
                ImLoglap4 = atan2(lap4->imag, lap4->real);

                gammaelement = ImLoglap1 + ImLoglap2 + ImLoglap3 + ImLoglap4;
                // Require gammaelement in [-PI, PI], since MATLAB imag(log(complex)) is in this range.
                while(gammaelement > PI){
                    gammaelement -= 2*PI;
                }
                while(gammaelement < -PI){
                    gammaelement += 2*PI;
                }

                gamma[m] += gammaelement;
            }
        }
        // Final Math
        gamma[m] /= 2*PI;
        gamma[m] += (double)1/Nphi;
    }

/************************************************************/

    // Turn gamma as CHERN_DOUBLE to gamma_int as CHERN_INT
    double gammasum = 0;
    int gammasum_int = 0;
    int gamma_int[Nphi];

    for(m = 0; m < N; m++){
        gammasum += gamma[m];
        gamma_int[m] = (int) floor(gamma[m]+0.5);//0.5 here for rounding
    }
    gammasum_int = (int) floor(gammasum+0.5); //0.5 here for rounding
//    printf("Sample Total Chern is %d\n", gammasum_int);

    // Report ErrChern Sum whihc != 1
    if(gammasum_int != 1)
        *howmanyerrchern += 1;

/*    // Print out CHERN_INT
    for(m = 0; m < N; m++){
        printf("%d\t", gamma_int[m]);
    }
    printf("\n");*/

/************************************************************/

    // Eigenvalue W
    double eigsum;
    double *eigvave;
    eigvave = (double *) malloc (Nphi * sizeof(double));

    double *eigv;
    for(i = 0; i < N; i++){
        eigsum = 0;
        for(counttheta1 = 0; counttheta1 < thetameshcut; counttheta1++){
            for(counttheta2 = 0; counttheta2 < thetameshcut; counttheta2++){
                eigv = eigvptrstore[counttheta1 * (thetameshcut + 1) + counttheta2];
                eigsum += eigv[i];
            }
        }
        eigvave[i] = (double) eigsum / (thetameshcut + 0) / (thetameshcut + 0); // Attention +0 here, because of the "<" sign.
    }

/************************************************************/

    FILE *output;

    // Record all states
    output = fopen("ensemble_data","a");
    for (i = 0; i < N; i++) {
        fprintf(output, "%7d%21.14f%6d\n", samplenum, eigvave[i], gamma_int[i]);
	//fprintf(output, "%7d%21.14f%10.4f\n", samplenum, eigvave[i], gamma[i]);
        if(gamma_int[i] != 0) *howmanyext += 1;
    }
    fclose(output);

/************************************************************/

    for(counttheta1 = 0; counttheta1 <= thetameshcut; counttheta1++){
        for(counttheta2 = 0; counttheta2 <= thetameshcut; counttheta2++){
            free(eigfptrstore[counttheta1 * (thetameshcut + 1) + counttheta2]);
            free(eigvptrstore[counttheta1 * (thetameshcut + 1) + counttheta2]);
        }
    }

/************************************************************/

    free(gamma);
    free(buffer_eigf00);
    free(buffer_eigf10);
    free(buffer_eigf01);
    free(buffer_eigf11);
    free(lap1);
    free(lap2);
    free(lap3);
    free(lap4);

    free(eigvave);
    return 0;
}
