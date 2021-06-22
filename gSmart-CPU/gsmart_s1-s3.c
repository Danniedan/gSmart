
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
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
	
  float milliseconds = 0;

	char * filename="/root/cyd/data/watdiv100_s1.txt";

  if ((f = fopen(filename, "r")) == NULL)
    exit(1);
	
	fscanf(f, "%d	%d	%ld", &N, &M, &nz);
	
  int num_row;
	num_row=N;
	int num_column;
	num_column=M;
	long nnz;
	nnz=nz;
	printf("num_row=%d, num_column=%d, nnz=%ld\n", num_row, num_column, nnz);
	
    /* reseve memory for matrices */

	int *subs=(int *) malloc(nnz * sizeof(int));
	int *obs=(int *) malloc(nnz * sizeof(int));
	int *pres=(int *) malloc(nnz*sizeof(int));
	
	long g;

  for (i=0; i<nnz; i++)
  {
		fscanf(f, "%d	%d	%d", &(subs[i]), &(pres[i]), &(obs[i]));
  }

	if (f !=stdin) fclose(f);
	
    /* store the matrix in csr format */
  
	int num_r;
	int *num_nonzeros=(int *) malloc(num_row * sizeof(int));
	for(i=0; i<num_row; i++)
		num_nonzeros[i]=0;
	int* Mr=(int *) malloc((num_row+1) * sizeof(int));

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

	int *Ap=(int *) malloc((num_r+1) * sizeof(int));
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
	
	int con = 7749058;
	int con_pre = 8;
	k=0;
	if( (Mr[con+1]-Mr[con])==1 )
	{
		//printf("nonempty row!\n");
		for(i=Ap[Mr[con]]; i<Ap[Mr[con]+1]; i++)
		{
			if( pres[i]==con_pre )
			{
				bind1[k] = obs[i];
				k++;
			}
		}
	}
	printf("k: %d\n", k);
	
	
	
	    /* evaluate 0th-level query edges */
  
  
	int *bind2 = (int *) malloc(8*k*maxl * sizeof(int));
	int edge_pre[8]={0, 1, 2, 6, 4, 3, 5, 7};
	int ind[8]={0, 0, 0, 0, 0, 0, 0, 0};
	int flag;
	int p=0;
	
	for(i=0; i<k; i++)
	{
		for(g=0; g<8; g++)
		{
			flag=0;
			if( (Mr[bind1[i]+1]-Mr[bind1[i]])==1 )
			{
				for(j=Ap[Mr[bind1[i]]]; j<Ap[Mr[bind1[i]]+1]; j++)
				{
					if( (pres[j]==edge_pre[g]) )
					{
						bind2[i*8*maxl+g*maxl+ind[g]] = obs[j];
						ind[g]++;
						flag=1;
					}
				}
				if(flag==0)
				{
					bind1[i]=-1;
					for(it=0; it<g; it++)
            ind[it] = 0;
					break;
				}
				
			}
		}
		if(flag==1)
    {
			for(g=0; g<8; g++)
				ind[g] = 0;
		}
		
	}

	return 0;
}

