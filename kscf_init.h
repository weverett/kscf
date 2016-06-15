#include "gaussian_shells.h"

void kscf_init(double *ex, double *cs, double *cp, double *cd, 
		 double *cf, double *cg, double *ch, double *ci, 
                 long int *kstart, long int *katom, long int *ktype, long int *kng, 
		 long int *kloc, long int *kmin, long int *kmax, 
	       long int *nshell, long int *nat, long int *num, long int *na, long int *nb, double *c, double *zaninp);


extern void kscf_init_(double *ex, double *cs, double *cp, double *cd, 
		 double *cf, double *cg, double *ch, double *ci, 
			 long int *kstart, long int *katom, long int *ktype, long int *kng, 
			 long int *kloc, long int *kmin, long int *kmax, 
  		       long int *nshell, long int *nat, long int *num, long int *na, long int *nb, double *c, double *zaninp)
{
  return kscf_init(ex,cs,cp,cd,cf,cg,ch,ci, 
                     kstart,katom,ktype,kng, 
		   kloc,kmin,kmax,nshell,nat,num,na,nb,c,zaninp);
}    

