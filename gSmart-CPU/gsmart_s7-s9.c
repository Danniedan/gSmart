
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>

const int num_con = 1;
const int num_pre = 3;


int main(int argc, char** argv)
{
  struct timeval startt;
  struct timeval endt;
  unsigned long timer;
	
  FILE *fr, *fc;
  int M;
	int N;
	long nnz, i, j, it, nzr, nzc;	
  
  float milliseconds = 0;

	char * filename1="./data/wat100r_s7.txt";
	char * filename2="./data/wat100c_s7.txt";

    if ((fr = fopen(filename1, "r")) == NULL) 
            exit(1);
	
	fscanf(fr, "%d	%d	%ld", &(M), &(N), &(nnz));

	int num_row;
	num_row=M;
	printf("num_row=%d, nnz=%ld\n", num_row, nnz);
	
    /* reseve memory for matrices */

	int *subs_r=(int *) malloc(nnz * sizeof(int));
	int* obs_r=(int *) malloc(nnz * sizeof(int));
	int* pres_r=(int *) malloc(nnz * sizeof(int));
	
     /* read raw data */
  
	int direc_consis[3]={0, 2, 3};
	int direc_op[1]={1};
	
	long g;
	int sub, pre, ob;
	g=0;

  for (i=0; i<nnz; i++)
  {
		fscanf(fr, "%d	%d	%d", &(sub), &(pre), &(ob));
		for(it=0; it<3; it++)
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
	printf("nzr=%d\n", nzr);

	if (fr !=stdin) fclose(fr);
	
  
    /* LSpM_CSR storing */
  
	int num_r;
	int *num_nonzeros = (int *) malloc(num_row * sizeof(int));
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
	printf("num_column=%d, nnz=%ld\n", num_column, nnz);
	
    /* reseve memory for matrices */

	int *obs_c = (int *) malloc(nnz * sizeof(int));
	int* subs_c = (int *) malloc(nnz * sizeof(int));
	int* pres_c = (int *) malloc(nnz * sizeof(int));
	
  
    /* read raw data */
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
	

    /* LSpM_CSC storing */
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
  
	int *bind1 = (int *) malloc(num_con*num_row * sizeof(int));
	for(i=0; i<num_con*num_row; i++)
		bind1[i] = 0;
	int *Bp1 = (int *) malloc(num_con * sizeof(int));
	for(i=0; i<(num_con); i++)
		Bp1[i]=0;
	
	
	//constants 27, 48, 
	int con[num_con] = {40};
	int con_pre[num_con] = {1};
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
	
	
	
	/* main compuation */
  	// evaluate 0th-level query edges

	int *bind2 = (int *) malloc(num_pre*k*maxl_r * sizeof(int));
	int edge_pre[num_pre]={0, 2, 3};
	int ind[num_pre];
	for(i=0; i<num_pre; i++)
		ind[i]=0;
	int p=0;
	
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

