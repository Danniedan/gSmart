
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>

const int num_con = 2;
const int num_pre = 1;
const int num_pre_op = 1;


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

	char * filename1="./data/wat100r_s10.txt";
	char * filename2="./data/wat100c_s10.txt";

    if ((fr = fopen(filename1, "r")) == NULL) 
            exit(1);
	
	fscanf(fr, "%d	%d	%ld", &(M), &(N), &(nnz));
	int num_row;
	num_row=M;
	printf("num_row=%d, nnz=%ld\n", num_row, nnz);
	
    /* reseve memory for matrices */

	int *subs_r = (int *) malloc(nnz * sizeof(int));
	int* obs_r = (int *) malloc(nnz * sizeof(int));
	int* pres_r = (int *) malloc(nnz * sizeof(int));
	
	
    /* read raw data for LSpM_CSR */
  
	int direc_consis[1]={2};
	int direc_op[3]={0, 1, 3};
	
	long g;
	int sub, pre, ob;
	g=0;
  
	for (i=0; i<nnz; i++)
	{
		fscanf(fr, "%d	%d	%d", &(sub), &(pre), &(ob));
		for(it=0; it<1; it++)
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
	
    /* reseve memory for matrices */

	int *obs_c = (int *) malloc(nnz * sizeof(int));
	int* subs_c = (int *) malloc(nnz * sizeof(int));
	int* pres_c = (int *) malloc(nnz * sizeof(int));
	

    /* read raw data for LSpM_CSC */
	g=0;
	for (i=0; i<nnz; i++)
	{
		fscanf(fc, "%d	%d	%d", &(sub), &(pre), &(ob));
		for(it=0; it<3; it++)
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
	printf("nzc=%d\n", nzc);
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
  
	int *bind1 = (int *) malloc(num_con*maxl_c * sizeof(int));
	for(i=0; i<num_con*maxl_c; i++)
		bind1[i] = 0;
	int *Bp1 = (int *) malloc(2*num_con * sizeof(int));
	for(i=0; i<(2*num_con); i++)
	  Bp1[i]=0;

	
	//constants 13387
	int con[num_con] = {13390, 13496};
	int con_pre[num_con] = {3, 1};
	int flag;
	if( ((Mc[con[0]+1]-Mc[con[0]])==1) && ((Mc[con[1]+1]-Mc[con[1]])==1) )
	{
		for(j=0; j<num_con; j++)
		{
			k=0;
			for(i=Ap_c[Mc[con[j]]]; i<Ap_c[Mc[con[j]]+1]; i++)
			{
				if( pres_c[i]==con_pre[j] )
				{
					bind1[j*maxl_c+k]=subs_c[i];
					k++;
				}
			}
			Bp1[j]=k;
			printf("Bp1[%d]=%d\n", j, Bp1[j]);
		}
	}
	k=Bp1[0];
	
	
	
	  /* main computation */
  
	int *bind2 = (int *) malloc(num_pre*k*maxl_r * sizeof(int));
	int edge_pre[num_pre]={2};
	int ind[num_pre];
	for(i=0; i<num_pre; i++)
		ind[i]=0;
	for(i=0; i<num_pre*k*maxl_r; i++)
		bind2[i]=-1;
	
	
	int *bind2_op = (int *) malloc( num_pre_op*k*maxl_c * sizeof(int));
	int edge_pre_op[num_pre_op]={0};
	int ind_op[num_pre_op];
	for(i=0; i<num_pre_op; i++)
		ind_op[i]=0;
	for(i=0; i<num_pre_op*k*maxl_r; i++)
		bind2_op[i]=-1;
	

    /* evaluate 0th-level query edges*/
  
	for(i=0; i<k; i++)
	{
		for(g=0; g<num_pre; g++)
		{
			flag=0;
			if( ((Mr[bind1[i]+1]-Mr[bind1[i]])==1)&&((Mc[bind1[i]+1]-Mc[bind1[i]])==1) )
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
			for(g=0; g<num_pre_op; g++)
			{
				flag=0;
				for(j=Ap_c[Mc[bind1[i]]]; j<Ap_c[Mc[bind1[i]]+1]; j++)
				{
					if( (pres_c[j]==edge_pre_op[g]) )
					{
						//printf("i=%d, pres_c[%d] = %d, edge_pre_op[%d] = %d\n", i, j, pres_c[j], g, edge_pre_op[g]);
						bind2_op[i*num_pre_op*maxl_c+g*maxl_c+ind_op[g]] = subs_c[j];
						ind_op[g]++;
						//if(ind_op[g]>1)
						//	printf("ind_op[%d]=%d\n", g, ind_op[g]);
						flag=1;
					}
				}
				if(flag==0)
				{
					bind1[i] = -1;
					for(it=0; it<g; it++)
						ind_op[it] = 0;
					break;
				}
			}
			
			if(flag==1)
			{
				for(g=0; g<num_pre; g++)
				{
					ind[g] = 0;
				}
				for(g=0; g<num_pre_op; g++)
				{
					ind_op[g] = 0;
				}
			}
		}
		
	}
	
    /* local tree-pruning */
  
	int *bind3 = (int *) malloc( 2*num_row * sizeof(int));
	for(i=0; i<2*num_row; i++)
		bind3[i] = 0;
		
	int p=0;
	
	for(i=0; i<k; i++)
	{
		if(bind1[i]!=-1)
		{
			bind3[bind1[i]]=1;
		}
	}
	for(i=0; i<Bp1[1]; i++)
	{
		bind3[num_row+bind1[maxl_c+i]]=1;
	}
	for(i=0; i<num_row; i++)
	{
		if((bind3[i]==1)&&(bind3[num_row+i]==1))
		{
			p++;
			
		}
	}
	printf("p=%d\n", p);
	
	

	return 0;
}

