#include "simint/simint_init.h"
#include "simint/eri/eri.h"

int *iscratch;
double *dscratch;
double *integrals;

void kscf_energx();

void oed__memory_ovl_batch_(int *NALPHA, int *NCOEF,
			    int *NCGTO1, int *NCGTO2,
			    int *NPGTO1, int *NPGTO2,
			    int *SHELL1, int *SHELL2,
			    double *X1, double *Y1, double *Z1,
			    double *X2, double *Y2, double *Z2,
			    double ALPHA[], double CC[],
			    int *SPHERIC,
			    
			    int *IMIN, int *IOPT,
			    int *ZMIN, int *ZOPT );

void oed__gener_ovl_batch_(int *IMAX, int *ZMAX,
			   int *NALPHA, int *NCOEF, int *NCSUM,
			   int *NCGTO1, int *NCGTO2,
			   int *NPGTO1, int *NPGTO2,
			   int *SHELL1, int *SHELL2,
			   double *X1, double *Y1, double *Z1,
			   double *X2, double *Y2, double *Z2,
			   double ALPHA[], double CC[], int CCBEG[], int CCEND[],
			   int *SPHERIC,
			   int *SCREEN,
			   int *ICORE,
			   
			   int *NBATCH,
			   int *NFIRST,
			   double *ZCORE);

void oed__memory_kin_batch_(int *NALPHA, int *NCOEF,
			    int *NCGTO1, int *NCGTO2,
			    int *NPGTO1, int *NPGTO2,
			    int *SHELL1, int *SHELL2,
			    double *X1, double *Y1, double *Z1,
			    double *X2, double *Y2, double *Z2,
			    double ALPHA[], double CC[],
			    int *SPHERIC,
			    
			    int *IMIN, int *IOPT,
			    int *ZMIN, int *ZOPT );

void oed__gener_kin_batch_(int *IMAX, int *ZMAX,
			   int *NALPHA, int *NCOEF, int *NCSUM,
			   int *NCGTO1, int *NCGTO2,
			   int *NPGTO1, int *NPGTO2,
			   int *SHELL1, int *SHELL2,
			   double *X1, double *Y1, double *Z1,
			   double *X2, double *Y2, double *Z2,
			   double ALPHA[], double CC[], int CCBEG[], int CCEND[],
			   int *SPHERIC,
			   int *SCREEN,
			   int *ICORE,
			   
			   int *NBATCH,
			   int *NFIRST,
			   double *ZCORE);

void oed__memory_nai_batch_(int *NALPHA, int *NCOEF,
			    int *NCGTO1, int *NCGTO2,
			    int *NPGTO1, int *NPGTO2,
			    int *SHELL1, int *SHELL2,
			    double *X1, double *Y1, double *Z1,
			    double *X2, double *Y2, double *Z2,
			    int *NUCLEI, double ALPHA[], double CC[],
			    int *SPHERIC,
			    
			    int *IMIN, int *IOPT,
			    int *ZMIN, int *ZOPT );

void oed__gener_nai_batch_(int *IMAX, int *ZMAX,
			   int *NALPHA, int *NCOEF, int *NCSUM,
			   int *NCGTO1, int *NCGTO2,
			   int *NPGTO1, int *NPGTO2,
			   int *SHELL1, int *SHELL2,
			   double *X1, double *Y1, double *Z1,
			   double *X2, double *Y2, double *Z2,
			   int *NUCLEI,
			   double XN[], double YN[], double ZN[], double NCHARGE[],
			   double ALPHA[], double CC[], int CCBEG[], int CCEND[],
			   int *SPHERIC,
			   int *SCREEN,
			   int *ICORE,
			   
			   int *NBATCH,
			   int *NFIRST,
			   double *ZCORE);


