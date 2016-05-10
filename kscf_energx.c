#include <stdio.h>
#include "kscf_energx.h"
#include "gaussian_shells.h"

void kscf_energx()
{
  
  
  //  printf("%f\n", firstSptr.x);

    printf("kwk inside kscf_energx\n");
    for(int i=0; i<10; ++i)
    {
      printf("s_shells[%d]\n", i);
      printf("x %f\n", s_shells[i].x);
      printf("y %f\n", s_shells[i].y);
      printf("z %f\n", s_shells[i].z);
      printf("\n");
    }

  for(int i=0; i<10; ++i)
    {
      printf("p_shells[%d]\n", i);
      printf("x %f\n", p_shells[i].x);
      printf("y %f\n", p_shells[i].y);
      printf("z %f\n", p_shells[i].z);
      printf("\n");
      }  
}

void kscf_energx_erd()
{
  printf("inside kscf_energx_erd\n");
  kscf_energx_oed();
}

void kscf_energx_opterd()
{
  printf("inside kscf_energx_opterd\n");
  kscf_energx_oed();
}

void kscf_energx_simint()
{
  printf("inside kscf_energx_simint\n");
  kscf_energx_oed();
}

void kscf_energx_oed()
{
  printf("inside kscf_energx_oed()\n");
}
