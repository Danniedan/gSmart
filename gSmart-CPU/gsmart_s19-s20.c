
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>


int main(int argc, char** argv)
{
  struct timeval startt;
  struct timeval endt;
  unsigned long timer;
	
  FILE *f, *fw;
  int M;
	int N;
	long nz, i, j, it;	
	
	char * filename="./data/wat100_s19.txt";

    if ((f = fopen(filename, "r")) == NULL) 
            exit(1);
	
	fscanf(f, "%d	%d	%ld", &N, &M, &nz);
	
  int num_row;
	num_row=N;
	int num_column;
	num_column=M;
	long nnz;
	nnz=nz;
	
    /* reseve memory for matrices */

	int *subs=(int *) malloc(nnz * sizeof(int));
	int *obs=(int *) malloc(nnz * sizeof(int));
	int *pres=(int *) malloc(nnz*sizeof(int));

	long g;

  
    /* read raw data */
  
  for (i=0; i<nnz; i++)
  {
    fscanf(f, "%d	%d	%d", &(subs[i]), &(pres[i]), &(obs[i]));
  }

	if (f !=stdin) fclose(f);
	
  
  
	  /* LSpM_CSR storing */
  
	int num_r;
	int *num_nonzeros = (int *) malloc(num_row * sizeof(int));
	for(i=0; i<num_row; i++)
		num_nonzeros[i]=0;
	int* Mr = (int *) malloc((num_row+1) * sizeof(int));
	
	int sum=0;
	int row_index=0;
	for(i=0;i<nnz;i++)
	{
		num_nonzeros[subs[i]]++;
	}
	free(subs);

	int count=0;
	int k;

	g=0;
	Mr[0]=0;
	for(i=0; i<num_row; i++)
	{
		if(num_nonzeros[i]==0)
			Mr[i+1]=g;
		else
		{	
			g++;
			Mr[i+1]=g;
		}
	}
	num_r=g;

  
	int* Ap = (int *) malloc((num_r+1) * sizeof(int));
	int maxl=0;
	Ap[0]=0;
	sum=0;
	g=0;

	for(i=0; i<num_row; i++)
	{
		if(num_nonzeros[i]!=0)
		{	
			g++;
			sum+=num_nonzeros[i];
			Ap[g]=sum;
			if(maxl < (Ap[g]-Ap[g-1]))
			{
				maxl=Ap[g]-Ap[g-1];	
			}
		}
	}
	free(num_nonzeros);
	
  
    
    /* light query evaluation */
  
	int *bind1 = (int *) malloc(maxl * sizeof(int));
  int *Bp1;
	
	//constants , 428086
	int con = 676407;
	int con_pre = 2;
	k=0;
	if( (Mr[con+1]-Mr[con])==1 )
	{
		for(i=Ap[Mr[con]]; i<Ap[Mr[con]+1]; i++)
		{
			if( pres[i]==con_pre )
			{
				bind1[k] = obs[i];
				k++;
			}
		}
	}
	


    /* main computation */
  
	int *bind2 = (int *) malloc(2*k*maxl * sizeof(int));
	int edge_pre[2]={0, 1};
	int ind[2]={0, 0};
	int flag;
	int p=0;
	
	//evaluate 0th-level query edges
	for(i=0; i<k; i++)
	{
		for(g=0; g<2; g++)
		{
			flag=0;
			if( (Mr[bind1[i]+1]-Mr[bind1[i]])==1 )
			{
				for(j=Ap[Mr[bind1[i]]]; j<Ap[Mr[bind1[i]]+1]; j++)
				{
					if( (pres[j]==edge_pre[g]) )
					{
						bind2[p*2*maxl+g*maxl+ind[g]] = obs[j];
						ind[g]++;
						flag=1;
					}
				}
				if(flag==0)
				{
          bind1[i] = -1;
          for(it=0; it<g; it++)
            ind[it] = 0;
          break;
				}
				
			}
		}
		if(flag==1)
		{
			p++;
			for(g=0; g<2; g++)
				ind[g] = 0;
		}
		
	}
  
	return 0;
}

