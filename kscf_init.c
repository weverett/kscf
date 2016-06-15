#include <stdio.h>
#include <stdlib.h> /*getenv, exit*/
#include <string.h> /*strcmp*/
#include "kscf_init.h"
#include "kscf_energx.h"
#include "gaussian_shells.h"
#include "simint/simint_init.h"
#include "simint/eri/eri.h"

int  totalNumShells;
int maxAm;

int primIndex[MAX_SHELLS];
int gaussIndex[MAX_SHELLS];

struct gaussian_shell *allShells[MAX_SHELLS];
struct gaussian_shell *shellPtr[MAX_AM];

struct atomInfo atomInfo;

struct gaussian_shell *sShells;
struct gaussian_shell *pShells;
struct gaussian_shell *dShells;
struct gaussian_shell *fShells;
struct gaussian_shell *gShells;
struct gaussian_shell *hShells;
struct gaussian_shell *iShells;

int numGauss;
int numOccAlpha;
int numOccBeta;
int numShells[MAX_AM];

void kscf_init(double *ex, double *cs, double *cp, double *cd, 
	       double *cf, double *cg, double *ch, double *ci, 
	       long int *kstart, long *katom, long int *ktype, long int *kng, 
	       long int *kloc, long int *kmin, long int *kmax, 
	       long int *nshell, long int *nat, long int *num, long int *na, long int *nb, double *c, double *zaninp)
{
  
  numGauss=*num;
  numOccAlpha=*na;
  numOccBeta=*nb;

  ///////////////////////////////////////////////////////////////////////////////
  // Grab KSCF_INTPKG environmental variable from rungms  
  // Must be ERD, OPTERD, or SIMINT  
  ///////////////////////////////////////////////////////////////////////////////
  
  char *envKSCF_INTPKG;
  
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

  for(int i=0; i<*nat; i++)
    {
      int coordStart=2*i+i;
      atomInfo.XN[i]=c[coordStart];
      atomInfo.YN[i]=c[coordStart+1];
      atomInfo.ZN[i]=c[coordStart+2];
      atomInfo.zCharge[i]=zaninp[i];
    }

  for(int i=0; i<*nshell; ++i)
    {
      if(ktype[i]==1)
	{
	  numShells[0]+=1;
	}
      else if(ktype[i]==2)
	{
	  if(kmin[i]==2)
	    {
	      numShells[1]+=1;
	    }
	  else if(kmin[i]==1)
	    {
	      numShells[0]+=1;
	      numShells[1]+=1;
	    }
	}
      else if(ktype[i]==3)
	{
	  numShells[2]+=1;
	}
      else if(ktype[i]==4)
	{
	  numShells[3]+=1;
	}
      else if(ktype[i]==5)
	{
	  numShells[4]+=1;
	}
      else if(ktype[i]==6)
	{
	  numShells[5]+=1;
	}
      else if(ktype[i]==7)
	{
	  numShells[6]+=1;
	}
    }
  maxAm=0;
  while(numShells[maxAm]!=0)
    {
      totalNumShells+=numShells[maxAm];
      maxAm+=1;
    }
  
  struct gaussian_shell s_shells[numShells[0]];
  struct gaussian_shell p_shells[numShells[1]];
  struct gaussian_shell d_shells[numShells[2]];
  struct gaussian_shell f_shells[numShells[3]];
  struct gaussian_shell g_shells[numShells[4]];
  struct gaussian_shell h_shells[numShells[5]];
  struct gaussian_shell i_shells[numShells[6]];

  shellPtr[0]=&s_shells[0];
  shellPtr[1]=&p_shells[0];
  shellPtr[2]=&d_shells[0];
  shellPtr[3]=&f_shells[0];
  shellPtr[4]=&g_shells[0];
  shellPtr[5]=&h_shells[0];
  shellPtr[6]=&i_shells[0];  


  
  //pointers to struct of all parameters for a given angular momentum  



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
  int gaussCounter=0;

  for(int i=0; i<*nshell; ++i)
    {
      int coordStart=3*(katom[i]-1);
      int shellCounter=sCounter+pCounter+dCounter+fCounter+gCounter+hCounter+iCounter;
      
      if(ktype[i]==1)
	{
	  allocate_gaussian_shell(kng[i], &s_shells[sCounter]);
	  allShells[shellCounter]=&s_shells[sCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[sCounter]=gaussCounter;
	  gaussCounter+=1;
	  
	  s_shells[sCounter].am=0;
	  s_shells[sCounter].nprim=kng[i];
	  s_shells[sCounter].x=c[coordStart];
	  s_shells[sCounter].y=c[coordStart+1];
	  s_shells[sCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      s_shells[sCounter].alpha[j]=ex[kstart[i]+j-1];
      	      s_shells[sCounter].coef[j]=cs[kstart[i]+j-1];
	      }

	  sCounter+=1;
	}
      else if(ktype[i]==2)
	{
	  if(kmin[i]==2)
  //ktype=kmin=2, so this is a standard P shell
	    {
	  allocate_gaussian_shell(kng[i], &p_shells[pCounter]);
	  allShells[shellCounter]=&p_shells[pCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[numShells[0]+3*pCounter]=gaussCounter;	  
	  gaussCounter+=3;
	  
	  p_shells[pCounter].am=1;
	  p_shells[pCounter].nprim=kng[i];
	  p_shells[pCounter].x=c[coordStart];
	  p_shells[pCounter].y=c[coordStart+1];
	  p_shells[pCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      p_shells[pCounter].alpha[j]=ex[kstart[i]+j-1];
	      p_shells[pCounter].coef[j]=cp[kstart[i]+j-1];
	      }

	  pCounter+=1;
	    }
	  else if(kmin[i]==1)
	    {
  //ktype=2 and kmin=1, so this is an L shell
	  allocate_gaussian_shell(kng[i], &s_shells[sCounter]);
	  allShells[shellCounter]=&s_shells[sCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[sCounter]=gaussCounter;	  	  
	  gaussCounter+=1;	  
	  allocate_gaussian_shell(kng[i], &p_shells[pCounter]);
  //shellCounter+1 to account for the just-added s_l component
	  allShells[shellCounter+1]=&p_shells[pCounter];
	  gaussIndex[shellCounter+1]=gaussCounter;
	  primIndex[numShells[0]+pCounter]=gaussCounter;	  	  	  
	  gaussCounter+=3;	  

	  s_shells[sCounter].am=0;
	  s_shells[sCounter].nprim=kng[i];
	  s_shells[sCounter].x=c[coordStart];
	  s_shells[sCounter].y=c[coordStart+1];
	  s_shells[sCounter].z=c[coordStart+2];	  
	  p_shells[pCounter].am=1;
	  p_shells[pCounter].nprim=kng[i];
	  p_shells[pCounter].x=c[coordStart];
	  p_shells[pCounter].y=c[coordStart+1];
	  p_shells[pCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      s_shells[sCounter].alpha[j]=ex[kstart[i]+j-1];
	      s_shells[sCounter].coef[j]=cs[kstart[i]+j-1];	      
	      p_shells[pCounter].alpha[j]=ex[kstart[i]+j-1];
	      p_shells[pCounter].coef[j]=cp[kstart[i]+j-1];
	      }	  

	  sCounter+=1;
	  pCounter+=1;
	    }
	}
      else if(ktype[i]==3)
	{
	  allocate_gaussian_shell(kng[i], &d_shells[dCounter]);
	  allShells[shellCounter]=&d_shells[dCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[numShells[0]+numShells[1]+dCounter]=gaussCounter;	  	  	  	  
	  gaussCounter+=6;	  

	  d_shells[dCounter].am=2;
	  d_shells[dCounter].nprim=kng[i];
	  d_shells[dCounter].x=c[coordStart];
	  d_shells[dCounter].y=c[coordStart+1];
	  d_shells[dCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      d_shells[dCounter].alpha[j]=ex[kstart[i]+j-1];
	      d_shells[dCounter].coef[j]=cd[kstart[i]+j-1];
	    }	  
	  
	  dCounter+=1;
	}
      else if(ktype[i]==4)
	{
	  allocate_gaussian_shell(kng[i], &f_shells[fCounter]);
	  allShells[shellCounter]=&f_shells[fCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[numShells[0]+numShells[1]+numShells[2]+fCounter]=gaussCounter;	  	  	  	  	  
	  gaussCounter+=7;	  

	  f_shells[fCounter].am=3;
	  f_shells[fCounter].nprim=kng[i];
	  f_shells[fCounter].x=c[coordStart];
	  f_shells[fCounter].y=c[coordStart+1];
	  f_shells[fCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      f_shells[fCounter].alpha[j]=ex[kstart[i]+j-1];
	      f_shells[fCounter].coef[j]=cf[kstart[i]+j-1];
	    }	  	  
	  
	  fCounter+=1;
	}
      else if(ktype[i]==5)
	{
	  allocate_gaussian_shell(kng[i], &g_shells[gCounter]);
	  allShells[shellCounter]=&g_shells[gCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[numShells[0]+numShells[1]+numShells[2]+numShells[3]+gCounter]=gaussCounter;	  	  	  	  	  	  
	  gaussCounter+=9;	  

	  g_shells[gCounter].am=4;
	  g_shells[gCounter].nprim=kng[i];
	  g_shells[gCounter].x=c[coordStart];
	  g_shells[gCounter].y=c[coordStart+1];
	  g_shells[gCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      g_shells[gCounter].alpha[j]=ex[kstart[i]+j-1];
	      g_shells[gCounter].coef[j]=cg[kstart[i]+j-1];
	    }	  	  	  
	  
	  gCounter+=1;
	}
      else if(ktype[i]==6)
	{
	  allocate_gaussian_shell(kng[i], &h_shells[hCounter]);
	  allShells[shellCounter]=&h_shells[hCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[numShells[0]+numShells[1]+numShells[2]+numShells[3]+numShells[4]+hCounter]=gaussCounter;	  	  	  	  	  	  	  
	  gaussCounter+=11;	  

	  h_shells[hCounter].am=5;
	  h_shells[hCounter].nprim=kng[i];
	  h_shells[hCounter].x=c[coordStart];
	  h_shells[hCounter].y=c[coordStart+1];
	  h_shells[hCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      h_shells[hCounter].alpha[j]=ex[kstart[i]+j-1];
	      h_shells[hCounter].coef[j]=ch[kstart[i]+j-1];
	    }	  	  	  
	  
	  hCounter+=1;
	}
      else if(ktype[i]==7)
	{
	  allocate_gaussian_shell(kng[i], &i_shells[iCounter]);
	  allShells[shellCounter]=&i_shells[iCounter];
	  gaussIndex[shellCounter]=gaussCounter;
	  primIndex[numShells[0]+numShells[1]+numShells[2]+numShells[3]+numShells[4]+numShells[5]+iCounter]=gaussCounter;	  	  	
	  gaussCounter+=13;	  

	  i_shells[iCounter].am=6;
	  i_shells[iCounter].nprim=kng[i];
	  i_shells[iCounter].x=c[coordStart];
	  i_shells[iCounter].y=c[coordStart+1];
	  i_shells[iCounter].z=c[coordStart+2];

	  for(int j=0; j<kng[i]; ++j)
	    {
	      i_shells[iCounter].alpha[j]=ex[kstart[i]+j-1];
	      i_shells[iCounter].coef[j]=ci[kstart[i]+j-1];
	    }	  	  	  	  
	  
	  iCounter+=1;
	}
    }
 struct gaussian_shell *sShells=&s_shells[0];
 struct gaussian_shell *pShells=&p_shells[0];
 struct gaussian_shell *dShells=&d_shells[0];
 struct gaussian_shell *fShells=&f_shells[0];
 struct gaussian_shell *gShells=&g_shells[0];
 struct gaussian_shell *hShells=&h_shells[0];
 struct gaussian_shell *iShells=&i_shells[0];  
  ///////////////////////////////////////////////////////////////////////////////
  //  Begin SCF energy computation
  ///////////////////////////////////////////////////////////////////////////////
  kscf_energx();
}












