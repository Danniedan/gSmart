
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

const int num_con = 1;
const int num_pre = 2;


__global__ void SPMV(int NUM_THREAD, int k, int maxl_r, int* Ap_r, int* Mr, int* pres_r, int* obs_r, int* bind1, int* bind2)
//__global__ void SPMV(int NUM_THREAD, int p, int k, int maxl, int* Ap, int* Mr, int* pres, int* obs, int* Bp2, int* bind1, int* bind2)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int edge_pre[num_pre]={1, 2};
	int ind[num_pre]={0, 0};
	if (idx < NUM_THREAD) {
		int i, j, g, it, flag;
		int start, end;
		start = k*idx/NUM_THREAD;
		end = k*(idx+1)/NUM_THREAD;
		//start = p+k*idx/NUM_THREAD;
		//end = p+k*(idx+1)/NUM_THREAD;
		
		for(i=start; i<end; i++)
		{
			//printf("i=%d\n", i);
		
			for(g=0; g<num_pre; g++)
			{
				flag=0;
				if( (Mr[bind1[i]+1]-Mr[bind1[i]])==1 )
				{
					//printf("sub=%d\n", bind[i]);
					for(j=Ap_r[Mr[bind1[i]]]; j<Ap_r[Mr[bind1[i]]+1]; j++)
					{
						if( (pres_r[j]==edge_pre[g]) )
						{
							//printf("pres_r[%d] = %d, edge_pre[%d] = %d\n", j, pres_r[j], g, edge_pre[g]);
							bind2[i*num_pre*maxl_r+g*maxl_r+ind[g]] = obs_r[j];
							ind[g]++;
							flag=1;
						}
					}
					if(flag==0)
					{
						bind1[i] = -1;
						//Bp1[0]--;
						for(it=0; it<g; it++)
							ind[it] = 0;
						break;
					}
					
				}
			}
			if(flag==1)
			{
				//printf("i=%d: ", i);
				//for(g=0; g<num_pre; g++)
				{
					//Bp2[g*k+i] = ind[g];
					//printf("Bp2[%d]=%d\n", 2*2*p+2*g+1, Bp2[2*2*p+2*g+1]);
				}
				for(g=0; g<num_pre; g++)
					ind[g] = 0;
			}	
		}
	}
}

int com(double* A, double* B, int num) {
	int flag = 1;
	for (int i = 0; i < num; i++) {
		if (fabs(A[i] - B[i]) > 0.00001) {
			printf("%lf %lf %d\n", A[i], B[i], i);
			flag = 0;
		}
	}
	return flag;
}

int compare(int* A, int* B, int num) {
	for (int i = 0; i < num; i++) {
		if (A[i] != B[i])
			return 0;
	}
	return 1;
}

int main(int argc, char** argv)
{
  struct timeval startt;
  struct timeval endt;
  unsigned long timer;

  FILE *fr, *fc;
  int M;
	int N;
	long nnz, i, j, it, nzr, nzc;	
	int round;
	
	
	char * filename1="./data/wat100r_s7.txt";
	char * filename2="./data/wat100c_s7.txt";

  if ((fr = fopen(filename1, "r")) == NULL) 
    exit(1);
	
	fscanf(fr, "%d	%d	%ld", &(M), &(N), &(nnz));
	
	int num_row;
	num_row=M;
	
    /* reseve memory for matrices */

	int *subs_r = (int *) malloc(nnz * sizeof(int));
	int* obs_r = (int *) malloc(nnz * sizeof(int));
	int* pres_r = (int *) malloc(nnz * sizeof(int));
	
	
    /* read raw data for LSpM_CSR */
  
	int direc_consis[2]={1, 2};
	int direc_op[1]={0};
	
	long g;
	int sub, pre, ob;
	g=0;

  for (i=0; i<nnz; i++)
  {
		fscanf(fr, "%d	%d	%d", &(sub), &(pre), &(ob));
		for(it=0; it<2; it++)
		{
			if(pre==direc_consis[it])
			{
				subs_r[g]=sub;
				pres_r[g]=pre;
				obs_r[g]=ob;
				g++;
				break;
			}
		}
  }
	nzr=g;
	
	if (fr !=stdin) fclose(fr);
	
	
    /* LSpM_CSR storing */
  
	int num_r;
	int *num_nonzeros=(int *) malloc(num_row * sizeof(int));
	for(i=0; i<num_row; i++)
		num_nonzeros[i]=0;
	int sum=0;
	int* Mr = (int *) malloc((num_row+1) * sizeof(int));
	
	for(i=0;i<nzr;i++)
	{
		num_nonzeros[subs_r[i]]++;
	}
	free(subs_r);

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

	
	int* Ap_r;
	cudaHostAlloc((int**)&Ap_r, (num_r+1) * sizeof(int), cudaHostAllocDefault);
	int maxl_r=0;
	Ap_r[0]=0;
	sum=0;
	g=0;

	for(i=0; i<num_row; i++)
	{
		if(num_nonzeros[i]!=0)
		{	
			g++;
			sum+=num_nonzeros[i];
			Ap_r[g]=sum;
			if(maxl_r < (Ap_r[g]-Ap_r[g-1]))
			{
				maxl_r=Ap_r[g]-Ap_r[g-1];	
			}
		}
	}
	free(num_nonzeros);
	
	
	
	if ((fc = fopen(filename2, "r")) == NULL) 
            exit(1);
	
	fscanf(fc, "%d	%d	%ld", &(M), &(N), &(nnz));
	
	int num_column;
	num_column=N;
	
    /* reseve memory for matrices */

	int *obs_c = (int *) malloc(nnz * sizeof(int));
	int* subs_c = (int *) malloc(nnz * sizeof(int));
	int* pres_c = (int *) malloc(nnz * sizeof(int));
  
  
    /* LSpM_CSC storing */
  
  g=0;
  for (i=0; i<nnz; i++)
  {
		fscanf(fc, "%d	%d	%d", &(sub), &(pre), &(ob));
		for(it=0; it<1; it++)
		{
			if(pre==direc_op[it])
			{
				subs_c[g]=sub;
				pres_c[g]=pre;
				obs_c[g]=ob;
				g++;
				break;
			}
		}
  }
	nzc=g;
	if (fc !=stdin) fclose(fc);
	
  
	int num_c;
	num_nonzeros = (int *) malloc(num_column * sizeof(int));
	for(i=0; i<num_column; i++)
		num_nonzeros[i]=0;
	int* Mc = (int *) malloc((num_column+1) * sizeof(int));
	
	
	for(i=0;i<nzc;i++)
	{
		num_nonzeros[obs_c[i]]++;
	}
	free(obs_c);


	g=0;
	Mc[0]=0;
	for(i=0; i<num_column; i++)
	{
		if(num_nonzeros[i]==0)
			Mc[i+1]=g;
		else
		{	
			g++;
			Mc[i+1]=g;
		}
	}
	num_c=g;

  
	int* Ap_c = (int *) malloc( (num_c+1) * sizeof(int));
	int maxl_c=0;
	Ap_c[0]=0;
	sum=0;
	g=0;

	for(i=0; i<num_column; i++)
	{
		if(num_nonzeros[i]!=0)
		{	
			g++;
			sum+=num_nonzeros[i];
			Ap_c[g]=sum;
			if(maxl_c < (Ap_c[g]-Ap_c[g-1]))
			{
				maxl_c=Ap_c[g]-Ap_c[g-1];	
			}
		}
	}
	free(num_nonzeros);
	
	
	
	  /*main computation */
  
	int *bind1 = (int *) malloc( num_con*num_row * sizeof(int));
	for(i=0; i<num_con*num_row; i++)
		bind1[i] = 0;
	int *Bp1 = (int *) malloc( num_con * sizeof(int));
	for(i=0; i<(num_con); i++)
		Bp1[i]=0;
	
	//contants 62, , 326
	int con[num_con] = {1265};
	int con_pre[num_con] = {0};
	int flag;
	k=0;
	if( ((Mc[con[0]+1]-Mc[con[0]])==1) )
	{
		for(j=0; j<num_con; j++)
		{
			k=0;
			for(i=Ap_c[Mc[con[j]]]; i<Ap_c[Mc[con[j]]+1]; i++)
			{
				if( pres_c[i]==con_pre[j] )
				{
					bind1[k] = subs_c[i];
					k++;
				}
			}
		}
	}
	printf("k: %d\n", k);
	
	
	
	int *bind2 = (int *) malloc( num_pre*k*maxl_r * sizeof(int));
	int edge_pre[num_pre]={1, 2};
	int ind[num_pre];
	for(i=0; i<num_pre; i++)
		ind[i]=0;
	int p=0;
	

	//evaluate 0th-level query edges
	for(i=0; i<k; i++)
	{
		for(g=0; g<num_pre; g++)
		{
			flag=0;
			if( (Mr[bind1[i]+1]-Mr[bind1[i]])==1 )
			{
				for(j=Ap_r[Mr[bind1[i]]]; j<Ap_r[Mr[bind1[i]]+1]; j++)
				{
					if( (pres_r[j]==edge_pre[g]) )
					{
						bind2[i*num_pre*maxl_r+g*maxl_r+ind[g]] = obs_r[j];
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
			for(g=0; g<num_pre; g++)
				ind[g] = 0;
		}
  }

	
	
	return 0;
}

