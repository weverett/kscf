#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "kscf_energx.h"
#include "gaussian_shells.h"
#include "mkl.h"

#define INDEX(i,j) (i>j) ? (i*numGauss+j) : (j*numGauss+i)

void kscf_energx()
{
  for(int i=0;i<10; i++)
    printf("i=%d numShells[i]=%d primIndex[i]=%d\n", i, numShells[i], primIndex[i]);

  int nnSize=numGauss*numGauss;
  double hCore[nnSize];
  double ovl_integrals[nnSize];
  double *ovl_integrals_batch;  
  //initialize integrals to zero, and overlap diagonals to 1
  int diagCounter=0;
  for(int i=0;i<nnSize;i++)
    {
	if(i==diagCounter)
	  {
	    ovl_integrals[i]=1;
	    diagCounter+=(numGauss+1);
	    
	  }
	else
	  {
	    ovl_integrals[i]=0;
	  }
	
      hCore[i]=0;
    }
  
  
  for(int i=0;i<totalNumShells;i++)
    {
      printf("i=%d gaussStart=%d\n", i, gaussIndex[i]);
    }

  for(int i=0;i<totalNumShells;i++)
    for(int j=i;j<totalNumShells;j++)
    {
      {
	printf("i=%d j=%d\n", i, j);
	//	NCOEFF=NALPHA and NCGTOx=NPGTOx as long as batching isn't implemented
	int NALPHA=allShells[i]->nprim + allShells[j]->nprim;
	int NCGTO1=1;
	int NCGTO2=1;
	int IMIN=0;
	int IOPT=0;
	int ZMIN=0;
	int ZOPT=0;
	int SPHERIC=0;
	int SCREEN=1;
	int buffer_offset=0;
	int ovl_nbatch=0;
	int kin_nbatch=0;
	int nai_nbatch=0;	
	int NUCLEI=3;
	int NCSUM=NCGTO1+NCGTO2;
	int ccbeg[2]={1, 1};
	int ccend[2]={allShells[i]->nprim, allShells[j]->nprim};

	//gaussSum=number of computed elements before this index
	int gaussSum=gaussIndex[i]*(gaussIndex[i]+1)/2;
	//batchGaussSum is the number of computed elements before an index in a single batch
	int batchGaussSum=0;	    		
	printf("gaussSum=%d\n",gaussSum);

	//compoundIndex is the index of the first element of this integral batch in the 1D array
	int compoundIndex=(gaussIndex[i]*numGauss)+gaussIndex[j];

	printf("compoundIndex=%d\n", compoundIndex);	

	double *kin_integrals_batch;
	double *nai_integrals_batch;

	double ALPHA[NALPHA];
	double CC[NALPHA];
	for(int k=0; k<allShells[i]->nprim; ++k)
	  {
	    ALPHA[k]=allShells[i]->alpha[k];
	    CC[k]=allShells[i]->coef[k];
	  }
	for(int k=allShells[i]->nprim; k<NALPHA; ++k)
	  {
	    ALPHA[k]=allShells[j]->alpha[k-allShells[i]->nprim];
	    CC[k]=allShells[j]->coef[k-allShells[i]->nprim];
	  }

	//total number of contractions is always 2 w/o batching

	
	printf("\nNALPHA=%d NPRIM(I)=%d NPRIM(J)=%d am(i)=%d am(j)=%d\n", NALPHA,allShells[i]->nprim, allShells[j]->nprim, allShells[i]->am, allShells[j]->am);
	//ovl mem

	
	if(i!=j)
	  {
	oed__memory_ovl_batch_(&NALPHA, &NALPHA,
			       &NCGTO1, &NCGTO2,			       
			       &allShells[i]->nprim, &allShells[j]->nprim,
			       &allShells[i]->am, &allShells[j]->am,
			       &allShells[i]->x, &allShells[i]->y, &allShells[i]->z,
			       &allShells[j]->x, &allShells[j]->y, &allShells[j]->z,			       			       			       
			       ALPHA, CC, &SPHERIC,

			       &IMIN, &IOPT,
			       &ZMIN, &ZOPT);

	iscratch = (int *)ALLOC(IOPT * sizeof(int));
	dscratch = (double *)ALLOC(ZOPT * sizeof(double));
	
	//ovl calc

	oed__gener_ovl_batch_(&IOPT, &ZOPT,
			      &NALPHA, &NALPHA, &NCSUM,
			      &NCGTO2, &NCGTO1,			       
			      &allShells[i]->nprim, &allShells[j]->nprim,
			      &allShells[i]->am, &allShells[j]->am,
			      &allShells[i]->x, &allShells[i]->y, &allShells[i]->z,
			      &allShells[j]->x, &allShells[j]->y, &allShells[j]->z,			       			       			       
			      ALPHA, CC, ccbeg, ccend,
			      &SPHERIC,
			      &SCREEN,
			      iscratch,

			      &ovl_nbatch,
			      &buffer_offset,
			      dscratch);

	ovl_integrals_batch=(double *)ALLOC(ovl_nbatch * sizeof(double));
	
	printf("ovl integrals, nbatch=%d\n", ovl_nbatch);
	if(ovl_nbatch>0)
	  {
	    for(int i=0;i<ovl_nbatch;i++)
	      {
		printf("dscratch[%d]=%f\n", i, dscratch[buffer_offset-1+i]);
		ovl_integrals_batch[i]=dscratch[buffer_offset-1+i];
	      }
	  }
	printf("\n");


	for(int k=0; k<2*(allShells[i]->am)+1; k++)
	  {
	    for(int l=0; l<2*(allShells[j]->am)+1; l++)
	      {
		if(k==0)
		  {
		    printf("in a row:%d\n", compoundIndex+l);
		    //elements in a row
		    ovl_integrals[compoundIndex+l]=ovl_integrals_batch[l];
		    printf("kwk10\n");
		    //		newIndex=compoundIndex+l+(k*numGauss-(gaussIndex[i]+k+1)/2+i);		
		  }
		else
		  {
		    //elements in a column
		    ovl_integrals[compoundIndex+k*numGauss]=ovl_integrals_batch[k];
		    //			newIndex=compoundIndex+l+batchGaussSum;
		  }
	      }
	  }	

/* i=0 gaussStart=0 */
/* i=1 gaussStart=1 */
/* i=2 gaussStart=2 */
/* i=3 gaussStart=5 */
/* i=4 gaussStart=6 */

		//		printf("newIndex=%d\n", newIndex);
		//		ovl_integrals[compoundIndex+i+(]=dscratch[buffer_offset-1+i];

	

	FREE(dscratch);
	FREE(iscratch);
	  }



	//kin mem
	oed__memory_kin_batch_(&NALPHA, &NALPHA,
			       &NCGTO1, &NCGTO2,			       
			       &allShells[i]->nprim, &allShells[j]->nprim,
			       &allShells[i]->am, &allShells[j]->am,
			       &allShells[i]->x, &allShells[i]->y, &allShells[i]->z,
			       &allShells[j]->x, &allShells[j]->y, &allShells[j]->z,			       			       			       
			       ALPHA, CC, &SPHERIC,

			       &IMIN, &IOPT,
			       &ZMIN, &ZOPT);

        iscratch = (int *)ALLOC(IOPT * sizeof(int));
	dscratch = (double *)ALLOC(ZOPT * sizeof(double));


	//kin calc
	oed__gener_kin_batch_(&IOPT, &ZOPT,
			      &NALPHA, &NALPHA, &NCSUM,
			      &NCGTO1, &NCGTO2,			       
			      &allShells[i]->nprim, &allShells[j]->nprim,
			      &allShells[i]->am, &allShells[j]->am,
			      &allShells[i]->x, &allShells[i]->y, &allShells[i]->z,
			      &allShells[j]->x, &allShells[j]->y, &allShells[j]->z,			       			       			       
			      ALPHA, CC, ccbeg, ccend,
			      &SPHERIC,
			      &SCREEN,
			      iscratch,

			      &kin_nbatch,
			      &buffer_offset,
			      dscratch);

	kin_integrals_batch=(double *)ALLOC(kin_nbatch * sizeof(double));	

	printf("kin integrals, nbatch=%d\n", kin_nbatch);	
	if(kin_nbatch>0)
	  {
	    for(int i=0;i<kin_nbatch;i++)
	      {
		printf("dscratch[%d]=%f\n", i, dscratch[buffer_offset-1+i]);
		kin_integrals_batch[i]=dscratch[buffer_offset-1+i];		
	      }
	  }
	printf("\n");


	
	FREE(dscratch);
	FREE(iscratch);	

	//nai mem
	oed__memory_nai_batch_(&NALPHA, &NALPHA,
			       &NCGTO1, &NCGTO2,			       
			       &allShells[i]->nprim, &allShells[j]->nprim,
			       &allShells[i]->am, &allShells[j]->am,
			       &allShells[i]->x, &allShells[i]->y, &allShells[i]->z,
			       &allShells[j]->x, &allShells[j]->y, &allShells[j]->z,			       			       			       
			       &NUCLEI, ALPHA, CC, &SPHERIC,
			       
			       &IMIN, &IOPT,
			       &ZMIN, &ZOPT);

	
	iscratch = (int *)ALLOC(IOPT * sizeof(int));
	dscratch = (double *)ALLOC(ZOPT * sizeof(double));	

	//nai calc
	oed__gener_nai_batch_(&IOPT, &ZOPT,
			      &NALPHA, &NALPHA, &NCSUM,
			      &NCGTO1, &NCGTO2,
			      &allShells[i]->nprim, &allShells[j]->nprim,
			      &allShells[i]->am, &allShells[j]->am,
			      &allShells[i]->x, &allShells[i]->y, &allShells[i]->z,
			      &allShells[j]->x, &allShells[j]->y, &allShells[j]->z,
			      &NUCLEI,
			      atomInfo.XN, atomInfo.YN, atomInfo.ZN, atomInfo.zCharge,
			      ALPHA, CC, ccbeg, ccend,
			      &SPHERIC,
			      &SCREEN,
			      iscratch,

			      &nai_nbatch,
			      &buffer_offset,
			      dscratch);

	nai_integrals_batch=(double *)ALLOC(nai_nbatch * sizeof(double));	

	printf("nai integrals, nbatch=%d\n", nai_nbatch);		
	if(nai_nbatch>0)
	  {
	    for(int i=0;i<nai_nbatch;i++)
	      {
		printf("dscratch[%d]=%f\n", i, dscratch[buffer_offset-1+i]);
		nai_integrals_batch[i]=dscratch[buffer_offset-1+i];		
	      }
	  }
	printf("\n");
	
	int kGauss=2*(allShells[i]->am)+1;
	int lGauss=2*(allShells[j]->am)+1;	  
	for(int k=0; k<kGauss; k++)
	  {
	    for(int l=0; l<lGauss;l++)
	      {
		if(k==0)
		  {
		    printf("in a row:%d\n", compoundIndex+l);
		    //elements in a row
		    hCore[compoundIndex+l]=kin_integrals_batch[l]+nai_integrals_batch[l];
		    //		newIndex=compoundIndex+l+(k*numGauss-(gaussIndex[i]+k+1)/2+i);		
		  }
		else if(k>=l)
		  {
		    //elements in a column
		    printf("in a col:%d\n", compoundIndex+k*numGauss+l);
		    hCore[compoundIndex+k*numGauss+l]=kin_integrals_batch[l*kGauss+k]+nai_integrals_batch[l*kGauss+k];		    
		    //			newIndex=compoundIndex+l+batchGaussSum;
		  }
	      }
	  }		



		FREE(kin_integrals_batch);
		FREE(nai_integrals_batch);
		FREE(dscratch);
		FREE(iscratch);
	/* for(int i=0; i<nnSize; i++) */
	/*   { */
	/*     printf("i=%d, hcore=%f\n", i, hCore[i]); */
	/*   } */


      }
    }

  //copy upper diagonal of hcore to lower diagonal
	for(int i=0; i<numGauss; i++)
	  {
	    for(int j=i;j<numGauss;j++)
	      {
		if(i!=j)
		  {
		hCore[j*numGauss+i]=hCore[i*numGauss+j];
		  }
	      }
	  }
	for(int i=0; i<nnSize; i++)
	  {
	    printf("i=%d, hcore=%f\n", i, hCore[i]);
	  }	
	////////////////////////////////////////////////////////////////////////////
	//  diagonalize overlap matrix
	////////////////////////////////////////////////////////////////////////////

	//FEAST is Fortran routine, so we specify our matrix as lower triangular
	char *uplo = "L";
	
	//we are diagonalizing entire NxN matrix, lda=N
	MKL_INT N = (MKL_INT) numGauss;
	MKL_INT lda = N;

	//array of feast settings
	MKL_INT fpm[128];
	feastinit(fpm);
	fpm[0]=1;

	//realtive error on the trace
	double epsout;

	//number of refinement loops
	MKL_INT loop;

	//lower/upper bound of eigenvalue search interval [emin,emax]
	double Emin=-100;
	double Emax=100;

	//initial guess for subspace dimension to be used
	MKL_INT M0=N;

	//array of eigenvalues
	double *E=(double *) malloc (sizeof(double)*N);

	//array of eigenvectors
	double X[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  X[i]=0.0;

	//total number of eigenvalues found in the interval
	MKL_INT M;

	//residual
	double *res=(double *) malloc (sizeof(double) * N);

	//diagonalization error code
	MKL_INT info=0;

	dfeast_syev(
		    uplo, //IN: UPLO = U, upper triangular matrix
		    &N, //IN: size of problem
		    ovl_integrals, //IN: dense matrix A
		    &lda, //INT: leading dimension of matrix A
		    fpm, //IN/OUT: array of MKL routine options
		    &epsout, //OUT: relative error of the trace
		    &loop, //OUT: number of refinement loops executed
		    &Emin, //IN: lower bound of search interval
		    &Emax, //IN: upper bound of search interval
		    &M0, //IN: initial guess for subspace dimension to be used
		    E, //OUT: the first M entries of eigenvalues
		    X, //IN/OUT: the first M entries of eigenvectors
		    &M, //OUT: the total number of eigenvalues found in the interval [emin,emax]
		    res, //OUT: the first M components contain the relative residual vector
		    &info //OUT:error code
		    );

	printf("eigenvalues\n");	
	for(MKL_INT i=0; i<N; i++)
	  {
	    printf("i=%lld, eval=%f\n", i, E[i]);
	  }
	//take the inverse square root of eigenvalues, in place
	vdInvSqrt(numGauss, E, E);

	printf("eigenvalues invsqrt\n");	
	for(MKL_INT i=0; i<N; i++)
	  {
	    printf("i=%lld, eval=%f\n", i, E[i]);
	  }	
	
	printf("eigenvectors\n");
	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("i=%lld, eval=%f\n", i, X[i]);
	  }

	//store X[i] for transpose in multiplication
	double Y[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  Y[i]=X[i];	
	
	//multiply X by invsqrt(eigenvalues)
	for(MKL_INT i=0; i<numGauss; i++)
	  {
	    for(MKL_INT j=0; j<numGauss; j++)
	      {
		X[i*numGauss+j]*=E[i];
	      }
	  }

	const double alpha=1.0;
	const double beta=1.0;

	//allocate space for symmetric orthogonalization matrix S
	double S[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  S[i]=0.0;
	
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
		      numGauss, numGauss, numGauss, alpha, X, numGauss, Y, numGauss, beta, S, numGauss);
	
	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("i=%lld, X=%f, Y=%f\n", i, X[i], Y[i]);
	  }

	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("i=%lld, S=%f, HCORE=%f\n", i, S[i], hCore[i]);
	  }

	
	////////////////////////////////////////////////////////////////////////////
	//  build the initial guess density
	////////////////////////////////////////////////////////////////////////////

	//temporary S^T * H_core storage
	double Ftmp[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  Ftmp[i]=0.0;

	//initial guess fock matrix in orthonormal AO basis	
	double F0[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  F0[i]=0.0;		

	// S^T * H_core = Htmp
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
		      numGauss, numGauss, numGauss, alpha, S, numGauss, hCore, numGauss, beta, Ftmp, numGauss);

	// Htmp * S = F0
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		      numGauss, numGauss, numGauss, alpha, Ftmp, numGauss, S, numGauss, beta, F0, numGauss);

	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("i=%lld, F0=%f\n", i, F0[i]);
	  }

	//diagonalize Fock matrix

	dfeast_syev(
		    uplo, //IN: UPLO = U, upper triangular matrix
		    &N, //IN: size of problem
		    F0, //IN: dense matrix A
		    &lda, //INT: leading dimension of matrix A
		    fpm, //IN/OUT: array of MKL routine options
		    &epsout, //OUT: relative error of the trace
		    &loop, //OUT: number of refinement loops executed
		    &Emin, //IN: lower bound of search interval
		    &Emax, //IN: upper bound of search interval
		    &M0, //IN: initial guess for subspace dimension to be used
		    E, //OUT: the first M entries of eigenvalues
		    X, //IN/OUT: the first M entries of eigenvectors
		    &M, //OUT: the total number of eigenvalues found in the interval [emin,emax]
		    res, //OUT: the first M components contain the relative residual vector
		    &info //OUT:error code
		    );

	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("initial evectors=%f\n",X[i]);
	  }

	//MO coefficients in non-orthogonal AO basis
	double C0[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  C0[i]=0.0;			

	// S * C'0 = C0
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		      numGauss, numGauss, numGauss, alpha, S, numGauss, X, numGauss, beta, C0, numGauss);	

	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("i=%lld; initial mo coefficients=%f\n",i,C0[i]);
	  }

	//density matrix
	double D0[N*N];
	for (MKL_INT i=0; i<N*N; i++)
	  D0[i]=0.0;


	//form the density matrix. strong parallellization target
  for(int i=0;i<numGauss;i++)
    for(int j=i;j<numGauss;j++)
    {
      {
	for(int k=0;k<numOccAlpha;k++)
	  {
	    D0[i*numGauss+j]+=C0[k*numGauss+i]*C0[k*numGauss+j];
	    //fill in lower triangular
	    D0[j*numGauss+i]=D0[i*numGauss+j];
	  }
      }
    }
	for(MKL_INT i=0; i<N*N; i++)
	  {
	    printf("i=%lld; initial density=%f\n",i,D0[i]);
	  }

	
	////////////////////////////////////////////////////////////////////////////
	//  compute the initial SCF energy
	////////////////////////////////////////////////////////////////////////////

	double e0=0.0;

	for(int i=0; i<nnSize; i++)
	  {
	    e0+=D0[i]*(2*hCore[i]);
	  }
	printf("initial electronic energy=%.10lf\n", e0);
	//kwk note: energy is wrong compared to project3 by 0.000X

	
  ////////////////////////////////////////////////////////////////////////////
  // 2ei
  ////////////////////////////////////////////////////////////////////////////
  
  simint_init();
  //normalize all shells
  for(int i=0; i<maxAm; i++)
    {
  normalize_gaussian_shells(numShells[i], shellPtr[i]);
    }
  struct multishell_pair shellPairs[20];
  int shellAm[20];
  int pairCounter=0;
    for(int i=0;i<maxAm;i++)
      {
	for(int j=i;j<maxAm;j++)
	  {
	    shellPairs[pairCounter]=create_multishell_pair(numShells[i], shellPtr[i], numShells[j], shellPtr[j]);
	    shellAm[pairCounter]=(2*i+1)*(2*j+1);
	    pairCounter+=1;
	  }
      }

   /* struct multishell_pair ss_pair = create_multishell_pair(4, s_shells, 4, s_shells); */
   /* struct multishell_pair ps_pair = create_multishell_pair(1, p_shells, 4, s_shells); */
   /* struct multishell_pair pp_pair = create_multishell_pair(1, p_shells, 1, p_shells); */
    int ncomputed=0;
    int ntotal=0;
    //s=1 p=3 d=6 f=10...
    int numComponents[]={1,3,6,10};
    double *targets = simint_allocate_target(1789*1789*1789);    

    for(int x=0;x<pairCounter;x++)
      for(int y=0;y<x+1;y++)
	{
	  ncomputed=simint_compute_eri(shellPairs[x], shellPairs[y], targets);
	  printf("i.nshell=%d i.nprim=%d am1=%d am2=%d\n", shellPairs[x].nshell12, shellPairs[x].nprim, shellPairs[x].am1, shellPairs[x].am2);
	  printf("j.nshell=%d j.nprim=%d am1=%d am2=%d\n", shellPairs[y].nshell12, shellPairs[y].nprim, shellPairs[y].am1, shellPairs[y].am2);
	  ncomputed*=shellAm[x]*shellAm[y];
	  int iBound=numComponents[shellPairs[x].am1]*numShells[shellPairs[x].am1];
	  int jBound=numComponents[shellPairs[x].am2]*numShells[shellPairs[x].am2];
	  int kBound=numComponents[shellPairs[y].am1]*numShells[shellPairs[y].am1];
	  int lBound=numComponents[shellPairs[y].am2]*numShells[shellPairs[y].am2];
	  printf("bounds i/j/k/l=%d %d %d %d\n", iBound, jBound, kBound, lBound);
	  for(int i=0;i<iBound;i++)
	    for(int j=0;j<jBound;j++)
	      for(int k=0;k<kBound;k++)
		for(int l=0;l<lBound;l++)	    
		  {
		    if(i>=j && k>=l && (i*(i+1)/2+j)>=(k*(k+1)/2+l))
		      {
			int index=pow(iBound,3.0)*i+pow(jBound,2.0)*j+kBound*k+l;
			printf("index=%d\n", index);
			printf("i,j,k,l = %d,%d,%d,%d val=%f\n", i,j,k,l,targets[index]);
		      }		    
		  }
	  
	}

    /* ncomputed=simint_compute_eri(ss_pair, ss_pair, targets+ntotal); */
    /* printf("Computed %d contracted ssss integrals\n", ncomputed); */
    /* ntotal += ncomputed; */

    /* ncomputed = simint_compute_eri(ps_pair, ss_pair, targets+ntotal); */
    /* ncomputed *= 3; */
    /* printf("Computed %d contracted psss integrals\n", ncomputed); */
    /* ntotal += ncomputed; */

    /* ncomputed = simint_compute_eri(ps_pair, ps_pair, targets+ntotal); */
    /* ncomputed *= 9; */
    /* printf("Computed %d contracted psps integrals\n", ncomputed); */
    /* ntotal += ncomputed; */

    /* ncomputed = simint_compute_eri(pp_pair, ss_pair, targets+ntotal); */
    /* ncomputed *= 9; */
    /* printf("Computed %d contracted ppss integrals\n", ncomputed); */
    /* ntotal += ncomputed; */

    /* ncomputed = simint_compute_eri(pp_pair, ps_pair, targets+ntotal); */
    /* ncomputed *= 27; */
    /* printf("Computed %d contracted ppps integrals\n", ncomputed); */
    /* ntotal += ncomputed; */

    /* ncomputed = simint_compute_eri(pp_pair, pp_pair, targets+ntotal); */
    /* ncomputed *= 81; */
    /* printf("Computed %d contracted pppp integrals\n", ncomputed); */
    /* ntotal += ncomputed; */

    /* printf("++Computed %d contracted integrals\n", ntotal); */

    /*     for(int i = 0; i < ntotal; i++) */
    /*         printf(" %d  %12.8e\n", i, targets[i]); */

	/* for(int i=0; i<4; i++) */
	/* for(int j=0; j<4; j++) */
	/* for(int k=0; k<4; k++) */
	/* for(int l=0; l<4; l++) */
	/*   { */
	/*     if(i>=j && k>=l && i*j>=k*l) */
	/*       { */
	/* 	printf("i,j,k,l = %d,%d,%d,%d\n", i,j,k,l); */
	/*       } */
	/*   } */
    
    
    ////////////////////////////////////////////////////////////////////////////
    //  compute the new fock matrix
    ////////////////////////////////////////////////////////////////////////////	
    
    
    double F[N*N];
    for(int i=0; i<nnSize; i++)
      {
	F[i]=hCore[i];
      }

    /* for(int i=0; i<4; i++) */
    /*   for(int j=0; j<4;j++) */
    /*   { */
    /* 	F[0]+=D0[i*7+j]*(2.0 * targets[i*4+j]-targets[16*i+j]); */
    /*   } */
    /* printf("f0 is %f\n", F[0]); */
    
    /* for(int i=0; i<numOccAlpha; i++) */
    /*   for(int j=0; j<numOccAlpha; j++) */
    /* 	  for(int k=0; k<numOccAlpha; k++) */
    /* 	    for(int l=0; l<numOccAlpha; l++) */
    for(int i=0; i<numGauss; i++)
      for(int j=0; j<numGauss; j++)
    	  for(int k=0; k<numGauss; k++)
    	    for(int l=0; l<numGauss; l++)
	      {
		int ij=INDEX(i,j);
		int kl=INDEX(k,l);
		
		int ijkl=INDEX(ij,kl);
		int ik=INDEX(i,k);
		int jl=INDEX(j,l);
		int ikjl=INDEX(ik,jl);
    
		F[ij]+=D0[kl]*(2*targets[ijkl]-targets[ikjl]);
	      }
    printf("newFock\n");
    for(int i=0; i<nnSize; i++)
      printf("i=%d, val=%f", i, F[i]);
    ////////////////////////////////////////////////////////////////////////////
    //  clean up
    ////////////////////////////////////////////////////////////////////////////
    

    for(int i=0; i<totalNumShells; i++)
      {
	free_gaussian_shell(*allShells[i]);
      }
    /* free_gaussian_shell(s_shells[0]); */
    /* free_gaussian_shell(s_shells[1]); */
    /* free_gaussian_shell(s_shells[2]); */
    /* free_gaussian_shell(s_shells[3]); */
    /* free_gaussian_shell(p_shells[0]); */
    /* free_multishell_pair(ss_pair); */
    /* free_multishell_pair(ps_pair); */
    /* free_multishell_pair(pp_pair); */
    simint_free_target(targets);
    simint_finalize();
    
}
