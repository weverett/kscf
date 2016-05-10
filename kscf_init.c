#include <stdio.h>
#include <stdlib.h> /*getenv, exit*/
#include <string.h> /*strcmp*/
#include "kscf_init.h"
#include "kscf_energx.h"
#include "simint/simint_init.h"
#include "simint/eri/eri.h"
#include "simint/shell/shell.h"

void kscf_init(double* ex, double* cs, double* cp, double* cd, 
	       double* cf, double* cg, double* ch, double* ci, 
	       long int* kstart, long int* katom, long int* ktype, long int* kng, 
	       long int* kloc, long int* kmin, long int* kmax, 
	       long int* nshell, double* c)
{
  int numSshells=0;
  int numPshells=0;
  int numDshells=0;
  int numFshells=0;
  int numGshells=0;
  int numHshells=0;
  int numIshells=0;
  
  int numLshells=0;

  int coordStart;
  
  int indexLshells[*nshell];

  char* envKSCF_INTPKG;

  //  simint_init();

  
  ///////////////////////////////////////////////////////////////////////////////
  // Grab KSCF_INTPKG environmental variable from rungms  
  // Must be ERD, OPTERD, or SIMINT  
  ///////////////////////////////////////////////////////////////////////////////

  
  envKSCF_INTPKG=getenv("KSCF_INTPKG");
  if(envKSCF_INTPKG!=NULL)
    {
      if(strcmp(envKSCF_INTPKG,"ERD") 
      && strcmp(envKSCF_INTPKG,"OPTERD") 
      && strcmp(envKSCF_INTPKG,"SIMINT"))
	{
	  printf("KSCF_INTPKG is set to an incorrect value: %s\n", envKSCF_INTPKG);
	  exit(1);
	}
    }
  else
    {
      printf("KSCF_INTPKG is not set!!!\n");
      exit(1);
    }

  
  ///////////////////////////////////////////////////////////////////////////////
  //  Count number of each shell angular momentum, then create shell structures
  ///////////////////////////////////////////////////////////////////////////////

  
  for(int i=0; i<*nshell; ++i)
    {
      if(ktype[i]==1)
	{
	  numSshells+=1;
	}
      else if(ktype[i]==2)
	{
	  if(kmin[i]==2)
	    {
	      numPshells+=1;
	    }
	  else if(kmin[i]==1)
	    {
	      numLshells+=1;
	      indexLshells[i]=1;
	    }
	}
      else if(ktype[i]==3)
	{
	  numDshells+=1;
	}
      else if(ktype[i]==4)
	{
	  numFshells+=1;
	}
      else if(ktype[i]==5)
	{
	  numGshells+=1;
	}
      else if(ktype[i]==6)
	{
	  numHshells+=1;
	}
      else if(ktype[i]==7)
	{
	  numIshells+=1;
	}
    }

  struct gaussian_shell s_shells[numSshells+numLshells];
  struct gaussian_shell p_shells[numSshells+numLshells];
  struct gaussian_shell d_shells[numSshells];
  struct gaussian_shell f_shells[numSshells];
  struct gaussian_shell g_shells[numSshells];
  struct gaussian_shell h_shells[numSshells];
  struct gaussian_shell i_shells[numSshells];


  ///////////////////////////////////////////////////////////////////////////////
  //  Fill shell structures with shell info
  ///////////////////////////////////////////////////////////////////////////////


  int sCounter=0;
  int pCounter=0;
  int dCounter=0;
  int fCounter=0;
  int gCounter=0;
  int hCounter=0;
  int iCounter=0;

  int lCounter=0;
  for(int i=0; i<*nshell; ++i)
    {
      coordStart=3*(katom[i]-1);
      if(ktype[i]==1)
	{
	  allocate_gaussian_shell(kng[i], &s_shells[sCounter]);

	  s_shells[sCounter].am=0;
	  s_shells[sCounter].nprim=kng[i];
	  s_shells[sCounter].x=c[coordStart];
	  s_shells[sCounter].y=c[coordStart+1];
	  s_shells[sCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      s_shells[sCounter].alpha[j]=ex[kstart[i]+j];
	      s_shells[sCounter].coef[j]=cs[kstart[i]+j];
	      }
	  
	  sCounter+=1;
	}
      else if(ktype[i]==2)
	{
	  if(kmin[i]==2)
	    {
	  allocate_gaussian_shell(kng[i], &p_shells[pCounter]);

	  p_shells[pCounter].am=0;
	  p_shells[pCounter].nprim=kng[i];
	  p_shells[pCounter].x=c[coordStart];
	  p_shells[pCounter].y=c[coordStart+1];
	  p_shells[pCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      p_shells[pCounter].alpha[j]=ex[kstart[i]+j];
	      p_shells[pCounter].coef[j]=cp[kstart[i]+j];
	      }
	  
	  pCounter+=1;
	    }
	  else if(kmin[i]==1)
	    {
	  allocate_gaussian_shell(kng[i], &s_shells[sCounter]);
	  allocate_gaussian_shell(kng[i], &p_shells[pCounter]);

	  s_shells[sCounter].am=0;
	  s_shells[sCounter].nprim=kng[i];
	  s_shells[sCounter].x=c[coordStart];
	  s_shells[sCounter].y=c[coordStart+1];
	  s_shells[sCounter].z=c[coordStart+2];	  
	  p_shells[pCounter].am=0;
	  p_shells[pCounter].nprim=kng[i];
	  p_shells[pCounter].x=c[coordStart];
	  p_shells[pCounter].y=c[coordStart+1];
	  p_shells[pCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      s_shells[sCounter].alpha[j]=ex[kstart[i]+j];
	      s_shells[sCounter].coef[j]=cs[kstart[i]+j];	      
	      p_shells[pCounter].alpha[j]=ex[kstart[i]+j];
	      p_shells[pCounter].coef[j]=cp[kstart[i]+j];
	      }	  

	  sCounter+=1;
	  pCounter+=1;
	    }
	}
      else if(ktype[i]==3)
	{
	  allocate_gaussian_shell(kng[i], &d_shells[dCounter]);

	  d_shells[dCounter].am=0;
	  d_shells[dCounter].nprim=kng[i];
	  d_shells[dCounter].x=c[coordStart];
	  d_shells[dCounter].y=c[coordStart+1];
	  d_shells[dCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      d_shells[dCounter].alpha[j]=ex[kstart[i]+j];
	      d_shells[dCounter].coef[j]=cd[kstart[i]+j];
	    }	  
	  
	  dCounter+=1;
	}
      else if(ktype[i]==4)
	{
	  allocate_gaussian_shell(kng[i], &f_shells[fCounter]);

	  f_shells[fCounter].am=0;
	  f_shells[fCounter].nprim=kng[i];
	  f_shells[fCounter].x=c[coordStart];
	  f_shells[fCounter].y=c[coordStart+1];
	  f_shells[fCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      f_shells[fCounter].alpha[j]=ex[kstart[i]+j];
	      f_shells[fCounter].coef[j]=cf[kstart[i]+j];
	    }	  	  
	  
	  fCounter+=1;
	}
      else if(ktype[i]==5)
	{
	  allocate_gaussian_shell(kng[i], &g_shells[gCounter]);

	  g_shells[gCounter].am=0;
	  g_shells[gCounter].nprim=kng[i];
	  g_shells[gCounter].x=c[coordStart];
	  g_shells[gCounter].y=c[coordStart+1];
	  g_shells[gCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      g_shells[gCounter].alpha[j]=ex[kstart[i]+j];
	      g_shells[gCounter].coef[j]=cg[kstart[i]+j];
	    }	  	  	  
	  
	  gCounter+=1;
	}
      else if(ktype[i]==6)
	{
	  allocate_gaussian_shell(kng[i], &h_shells[hCounter]);

	  h_shells[hCounter].am=0;
	  h_shells[hCounter].nprim=kng[i];
	  h_shells[hCounter].x=c[coordStart];
	  h_shells[hCounter].y=c[coordStart+1];
	  h_shells[hCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      h_shells[hCounter].alpha[j]=ex[kstart[i]+j];
	      h_shells[hCounter].coef[j]=ch[kstart[i]+j];
	    }	  	  	  
	  
	  hCounter+=1;
	}
      else if(ktype[i]==7)
	{
	  allocate_gaussian_shell(kng[i], &i_shells[iCounter]);

	  i_shells[iCounter].am=0;
	  i_shells[iCounter].nprim=kng[i];
	  i_shells[iCounter].x=c[coordStart];
	  i_shells[iCounter].y=c[coordStart+1];
	  i_shells[iCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      i_shells[iCounter].alpha[j]=ex[kstart[i]+j];
	      i_shells[iCounter].coef[j]=ci[kstart[i]+j];
	    }	  	  	  	  
	  
	  iCounter+=1;
	}
    }

    for(int i=0; i<sCounter; ++i)
    {
      printf("s_shells[%d]\n", i);
      printf("x %f\n", s_shells[i].x);
      printf("y %f\n", s_shells[i].y);
      printf("z %f\n", s_shells[i].z);
      printf("\n");
    }

  for(int i=0; i<pCounter; ++i)
    {
      printf("p_shells[%d]\n", i);
      printf("x %f\n", p_shells[i].x);
      printf("y %f\n", p_shells[i].y);
      printf("z %f\n", p_shells[i].z);
      printf("\n");
      }
  /* struct gaussian_shell* firstSptr = &s_shells[0]; */
  /* struct gaussian_shell* firstPptr = &p_shells[0]; */
  /* struct gaussian_shell* firstDptr = &d_shells[0]; */
  /* struct gaussian_shell* firstFptr = &f_shells[0]; */
  /* struct gaussian_shell* firstGptr = &g_shells[0]; */
  /* struct gaussian_shell* firstHptr = &h_shells[0]; */
  /* struct gaussian_shell* firstIptr = &i_shells[0];   */
  kscf_energx();
}












