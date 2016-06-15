#ifndef GAUSSIAN_SHELLS_H
#define GAUSSIAN_SHELLS_H

#include "simint/simint_init.h"
#include "simint/eri/eri.h"

#define MAX_SHELLS 1000
#define MAX_ATOMS 1000
#define MAX_AM 7

/* defined in kscf_init.c */
/* extern int numSshells; */
/* extern int numPshells; */
/* extern int numDshells; */
/* extern int numFshells; */
/* extern int numGshells; */
/* extern int numHshells; */
/* extern int numIshells; */


extern struct gaussian_shell s_shells[];
extern struct gaussian_shell p_shells[];
extern struct gaussian_shell d_shells[];
extern struct gaussian_shell f_shells[];
extern struct gaussian_shell g_shells[];
extern struct gaussian_shell h_shells[];
extern struct gaussian_shell i_shells[];

extern struct gaussian_shell *sShells;
extern struct gaussian_shell *pShells;
extern struct gaussian_shell *dShells;
extern struct gaussian_shell *fShells;
extern struct gaussian_shell *gShells;
extern struct gaussian_shell *hShells;
extern struct gaussian_shell *iShells;

extern struct gaussian_shell *shellPtr[MAX_AM];
extern int totalNumShells;
extern int numGauss;
extern int numOccAlpha;
extern int numOccBeta;
extern int gaussIndex[MAX_SHELLS];
extern int primIndex[MAX_SHELLS];
extern int numShells[MAX_AM];
extern int maxAm;

extern struct gaussian_shell *allShells[MAX_SHELLS];

struct atomInfo
{
  double XN[MAX_ATOMS];
  double YN[MAX_ATOMS];
  double ZN[MAX_ATOMS];
  double zCharge[MAX_ATOMS];
};

extern struct atomInfo atomInfo;


#endif
