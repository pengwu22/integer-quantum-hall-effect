/********************************************************
    FUNCTION! UN-EXECUTABLE ALONE!

    Prepare for samplefunc.c

 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // Record Running Time

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


int presamplefunc(ptr_complex_double Coeff234, ptr_complex_double Coeff45, int *Coeff6, int Nphi, int quantaConcentration, double impurConcentration, int thetameshcut, int offhead);
int presamplefunc(ptr_complex_double Coeff234, ptr_complex_double Coeff45, int *Coeff6, int Nphi, int quantaConcentration, double impurConcentration, int thetameshcut, int offhead)
{
/*	int Nphi = 16;
	int quantaConcentration = 1;
	double impurConcentration = 10;
	int thetameshcut = 10;
	int offhead = 10;
*/
    double Area;
    Area = (double) Nphi / quantaConcentration;

    double L1;
    L1 = sqrt(Area);
    double L2 = L1;

	int m=0,n=0;
	int l=0,k=0;

    // Main Parameters: INPUT
    //int thetameshcut = 10; // !! Var4
    double thetastep;
    thetastep = 2 * PI / thetameshcut;
    //int offhead = 10; // Var5

    // Matrix Part: Size N
    int N = Nphi;

////////////////////////////////////////////////

    int dim_n, dim_nm, dim_ml, dim_nm1, dim_ml2;
    dim_n = offhead+1+offhead;
    dim_nm = dim_n * dim_n;
    dim_nm1 = dim_nm * (thetameshcut + 1);
    dim_ml = dim_n * Nphi;
    dim_ml2 = dim_ml * (thetameshcut + 1);

    ptr_complex_double buffer_complex;
    buffer_complex = (ptr_complex_double) malloc(sizeof(_complex_double));

    double buffer_Coeff2=0;
    double buffer_double=0;

    int counttheta1 = 0;
    int counttheta2 = 0;
    double theta1 = 0;
    double theta2 = 0;


    // Coeff6: n, k, l
    for(l = 1; l <= N; l++)
        for(k = l; k <= N; k++)
                    for(n = -offhead; n <= offhead; n++)
                        Coeff6[n+offhead+dim_n*( l+(k-1)*k/2-1 )] = math_LargeKroDeltaPBC(l,k-n,Nphi);// /L1/L2 to Coeff2

    // Coeff234: n, m, theta1
    theta1 = 0;
    for(counttheta1 = 0; counttheta1 <= thetameshcut; counttheta1++){

                    for(m = -offhead; m <= offhead; m++){
                        for(n = -offhead; n <= offhead; n++){
        // Prefactor II: Real double Exp
                            buffer_Coeff2 = exp(-PI*((double)m*m*L2/L1 + (double)n*n*L1/L2)/(2*Nphi))/L1/L2; // /L1/L2 from Coeff6, but really the normalization in Vmn
        // Prefactor III, IV-1: Complex Exp
                            buffer_double = (double)(-PI*m*n -n*theta1)/Nphi;
                            math_expimag(buffer_double, buffer_complex);// buffer_complex is re-valued here.
        // Prefactor C2*C34 // Attention to the initial value of buffer_complex
                            buffer_complex->real *= buffer_Coeff2;
                            buffer_complex->imag *= buffer_Coeff2;

                            Coeff234[(n+offhead) + dim_n*(m+offhead) + dim_nm*counttheta1] = *buffer_complex;
                        }
                    }
        ///////////////////////
        theta1 += thetastep;
    }

    // Coeff345: m, l, theta2
    theta2 = 0;
    for(counttheta2 = 0; counttheta2 <= thetameshcut; counttheta2++){

        for(l = 1; l <= N; l++){

            for(m = -offhead; m <= offhead; m++){
// Prefactor IV-2, V: Complex Exp
                buffer_double = (double)( m*theta2 - 2*PI*l*m )/Nphi;
                math_expimag(buffer_double, buffer_complex);// buffer_complex is re-valued here.

                Coeff45[(m+offhead) + dim_n*(l-1) + dim_ml*counttheta2] = *buffer_complex;
            }

        }
        ///////////////////////
        theta2 += thetastep;
    }

    free(buffer_complex);

	return 0;
}
