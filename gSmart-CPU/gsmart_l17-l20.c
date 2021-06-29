
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
const int num_pre_l2 = 1;
const int num_pre_l3 = 1;


__global__ void SPMV(int NUM_THREAD, int k, int maxl_r, int maxl_c, int num_row, int* Ap_r, int* Mr, int* Ap_c, int* Mc, 
	int* pres_r, int* obs_r, int* pres_c, int* subs_c, int* bind1, int* bind2, int* bind3)
{
	int i, j, g, it, flag, p;
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	
	int edge_pre_l2[num_pre_l2]={1};
	int ind_l2[num_pre_l2];
	for(i=0; i<num_pre_l2; i++)
		ind_l2[i]=0;
	
	int edge_pre_l3[num_pre_l3]={2};
	int ind_l3[num_pre_l3];
	for(i=0; i<num_pre_l3; i++)
		ind_l3[i]=0;
	
	if (idx < NUM_THREAD) {
		int start, end;
		start = k*idx/NUM_THREAD;
		end = k*(idx+1)/NUM_THREAD;
		//start = p+k*idx/NUM_THREAD;
		//end = p+k*(idx+1)/NUM_THREAD;
		
		for(i=start; i<end; i++)
		{
			//printf("i=%d\n", i);
			for(g=0; g<num_pre_l2; g++)
			{
				flag=0;
				if( (Mc[bind1[i]+1]-Mc[bind1[i]])==1 )
				{
					//printf("sub=%d\n", bind[i]);
					for(j=Ap_c[Mc[bind1[i]]]; j<Ap_c[Mc[bind1[i]]+1]; j++)
					{
						if( (pres_c[j]==edge_pre_l2[g]) )
						{
							//printf("pres_r[%d] = %d, edge_pre_l2[%d] = %d\n", j, pres_r[j], g, edge_pre_l2[g]);
							bind2[g*k*maxl_c+i*maxl_c+ind_l2[g]] = subs_c[j];
							ind_l2[g]++;
							flag=1;
						}
					}
					if(flag==0)
					{
						bind1[i] = -1;
						for(it=0; it<g; it++)
							ind_l2[it] = 0;
						break;
					}
					
				}
			}
			if(flag==1)
			{
				//计算1th-level query edges
				//for each result of v0 when v3=i
				//printf("i=%d, ind_l2[0]=%d\n", i, ind_l2[0]);
				for(p=0; p<ind_l2[0]; p++)
				{
					for(g=0; g<num_pre_l3; g++)
					{
						flag=0;
						if( (Mr[bind2[i*maxl_r+p]+1]-Mr[bind2[i*maxl_r+p]])==1 )
						{
							//printf("sub=%d\n", bind[i]);
							for(j=Ap_r[Mr[bind2[i*maxl_r+p]]]; j<Ap_r[Mr[bind2[i*maxl_r+p]]+1]; j++)
							{
								if( (pres_r[j]==edge_pre_l3[g]) )
								{
									//printf("pres_r[%d] = %d, edge_pre[%d] = %d\n", j, pres_r[j], g, edge_pre[g]);
									bind3[g*k*maxl_c*maxl_r+i*maxl_c*maxl_r+ind_l2[0]*maxl_r+ind_l3[g]] = obs_r[j];
									ind_l3[g]++;
									flag=1;
								}
							}
							if(flag==0)
							{
								bind2[i*maxl_c+p] = -1;
								for(it=0; it<g; it++)
									ind_l3[it] = 0;
								break;
							}
							
						}
					}
					if(flag==1)
					{
						for(it=0; it<g; it++)
							ind_l3[it] = 0;
					}
				}
				if(flag==0)
				{
					bind1[i] = -1;
					for(it=0; it<g; it++)
						ind_l2[it] = 0;
				}
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

	char * filename1="./data/wat100r_l17.txt";
	char * filename2="./data/wat100c_l17.txt";

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
	
	int direc_consis[2]={0, 2};
	int direc_op[1]={1};
	
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
	int *num_nonzeros = (int *) malloc(num_row * sizeof(int));
	for(i=0; i<num_row; i++)
		num_nonzeros[i]=0;
	int* Mr = (int *) malloc((num_row+1) * sizeof(int));
	
	int sum=0;
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
  
	int* Ap_r = (int *) malloc((num_r+1) * sizeof(int));
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
	

    /* read raw data for LSpM_CSC */
	
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
	
  
    /*LSpM_CSC storing */
  
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
  
	int* Ap_c = (int *) malloc((num_c+1) * sizeof(int));
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
	


	  /* light query evaluation */
	int *bind1 = (int *) malloc(maxl_r * sizeof(int));
	for(i=0; i<maxl_r; i++)
		bind1[i] = 0;
	
	//evaluation 184, 169, 52
	int con[num_con] = {16};
	int con_pre[num_con] = {0};
	int flag;
	if( (Mr[con[0]+1]-Mr[con[0]])==1 )
	{
		k=0;
		for(i=Ap_r[Mr[con[0]]]; i<Ap_r[Mr[con[0]]+1]; i++)
		{
			if( pres_r[i]==con_pre[0] )
			{
				bind1[k] = obs_r[i];
				k++;
			}
		}	
	}
	
	
  
    /* main computation */
  
	int *bind2 = (int *) malloc(num_pre_l2*k*maxl_c * sizeof(int));
	int edge_pre_l2[num_pre_l2]={1};
	int ind_l2[num_pre_l2];
	for(i=0; i<num_pre_l2; i++)
		ind_l2[i]=0;
	
	int *bind3 = (int *) malloc(num_pre_l3*k*maxl_c*maxl_r * sizeof(int));
	int edge_pre_l3[num_pre_l3]={2};
	int ind_l3[num_pre_l3];
	for(i=0; i<num_pre_l3; i++)
		ind_l3[i]=0;
		
	int p, q;
	
	//evaluate the 0th-level query edges
	for(i=0; i<k; i++)
	{
		for(g=0; g<num_pre_l2; g++)
		{
			flag=0;
			if( (Mc[bind1[i]+1]-Mc[bind1[i]])==1 )
			{
				for(j=Ap_c[Mc[bind1[i]]]; j<Ap_c[Mc[bind1[i]]+1]; j++)
				{
					if( (pres_c[j]==edge_pre_l2[g]) )
					{
						bind2[g*k*maxl_c+i*maxl_c+ind_l2[g]] = subs_c[j];
						ind_l2[g]++;
						flag=1;
					}
				}
				if(flag==0)
				{
					bind1[i] = -1;
					for(it=0; it<g; it++)
						ind_l2[it] = 0;
					break;
				}
				
			}
		}
		if(flag==1)
		{
			//evaluate the 1st-level query edges
			for(p=0; p<ind_l2[0]; p++)
			{
				for(g=0; g<num_pre_l3; g++)
				{
					flag=0;
					if( (Mr[bind2[i*maxl_r+p]+1]-Mr[bind2[i*maxl_r+p]])==1 )
					{
						for(j=Ap_r[Mr[bind2[i*maxl_r+p]]]; j<Ap_r[Mr[bind2[i*maxl_r+p]]+1]; j++)
						{
							if( (pres_r[j]==edge_pre_l3[g]) )
							{
								bind3[g*k*maxl_c*maxl_r+i*maxl_c*maxl_r+ind_l2[0]*maxl_r+ind_l3[g]] = obs_r[j];
								ind_l3[g]++;
								flag=1;
							}
						}
						if(flag==0)
						{
							bind2[i*maxl_c+p] = -1;
							for(it=0; it<g; it++)
								ind_l3[it] = 0;
							break;
						}
						
					}
				}
				if(flag==1)
				{
					for(it=0; it<g; it++)
						ind_l3[it] = 0;
				}
			}
			if(flag==0)
			{
				bind1[i] = -1;
				for(it=0; it<g; it++)
					ind_l2[it] = 0;
			}
		}
		
	}

	
	return 0;
}

