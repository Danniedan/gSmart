
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>

const int num_con = 2;
const int num_pre_l2 = 3;
const int num_pre_l3 = 1;


int main(int argc, char** argv)
{
  
  struct timeval startt;
  struct timeval endt;
  unsigned long timer;
  
  FILE *fr, *fc;
  int M;
	int N;
	long nnz, i, j, it, nzr, nzc;	
	
	char * filename1="./data/wat100r_f1.txt";
	char * filename2="./data/wat100c_f1.txt";

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
	
	int direc_consis[4]={0, 1, 3, 4};
	int direc_op[2]={1, 2};
	
	long g;
	int sub, pre, ob;
	g=0;
  
  for (i=0; i<nnz; i++)
  {
		fscanf(fr, "%d	%d	%d", &(sub), &(pre), &(ob));
		for(it=0; it<4; it++)
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
		for(it=0; it<2; it++)
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
	

	  /*light query evaluation */
  
	int *bind1_1 = (int *) malloc(maxl_c * sizeof(int));
	for(i=0; i<maxl_c; i++)
		bind1_1[i] = 0;
	int *bind1_2 = (int *) malloc(num_row * sizeof(int));
	for(i=0; i<num_row; i++)
		bind1_2[i] = 0;
	int *Bp1 = (int *) malloc(num_con * sizeof(int));
	for(i=0; i<(num_con); i++)
		Bp1[i]=0;
	
	//constants 2548749, 2549011, 2549079
	int con[num_con] = {107567, 107361};
	int con_pre[num_con] = {1, 2};
	int flag;
	if( (Mc[con[0]+1]-Mc[con[0]])==1 )
	{
		k=0;
		for(i=Ap_c[Mc[con[0]]]; i<Ap_c[Mc[con[0]]+1]; i++)
		{
			if( pres_c[i]==con_pre[0] )
			{
				bind1_1[k] = subs_c[i];
				k++;
			}
		}
		Bp1[0]=k;
		
	}
	
	if( (Mc[con[1]+1]-Mc[con[1]])==1 )
	{
		k=0;
		for(i=Ap_c[Mc[con[1]]]; i<Ap_c[Mc[con[1]]+1]; i++)
		{
			if( pres_c[i]==con_pre[1] )
			{
				bind1_2[subs_c[i]] = 1;
				k++;
			}
		}
		Bp1[1]=k;
		
	}
	
	
	
    /* main computation */
  
	int *bind2 = (int *) malloc(num_pre_l2*Bp1[0]*maxl_r * sizeof(int));
	int edge_pre_l2[num_pre_l2]={3, 4, 0};
	int ind_l2[num_pre_l2];
	for(i=0; i<num_pre_l2; i++)
		ind_l2[i]=0;
	for(i=0; i<num_pre_l2*Bp1[0]*maxl_r; i++)
		bind2[i]=-1;
	
	int *bind3 = (int *) malloc(num_pre_l3*Bp1[0]*maxl_r*maxl_r * sizeof(int));
	int edge_pre_l3[num_pre_l3]={1};
	int ind_l3[num_pre_l3];
	for(i=0; i<num_pre_l3; i++)
		ind_l3[i]=0;
	for(i=0; i<num_pre_l3*Bp1[0]*maxl_r*maxl_r; i++)
		bind3[i]=-1;
	
	int *tree1 = (int *) malloc(Bp1[0]*3*2 * sizeof(int));
	for(i=0; i<Bp1[0]*3*2; i++)
		tree1[i]=-1;
	int *tree2 = (int *) malloc(Bp1[0]*1*2 * sizeof(int));
	for(i=0; i<Bp1[0]*1*2; i++)
		tree2[i]=-1;
	int *tree3 = (int *) malloc(Bp1[0]*3*1 * sizeof(int));
	for(i=0; i<Bp1[0]*3*1; i++)
		tree3[i]=-1;
		
	int p, q;

	//evaluate the 0th-level query edges
	for(i=0; i<Bp1[0]; i++)
	{
		for(g=0; g<num_pre_l2; g++)
		{
			flag=0;
			if( (Mr[bind1_1[i]+1]-Mr[bind1_1[i]])==1 )
			{
				for(j=Ap_r[Mr[bind1_1[i]]]; j<Ap_r[Mr[bind1_1[i]]+1]; j++)
				{
					if( (pres_r[j]==edge_pre_l2[g]) )
					{
						bind2[g*Bp1[0]*maxl_r+i*maxl_r+ind_l2[g]] = obs_r[j];
						ind_l2[g]++;
						flag=1;
					}
				}
				if(flag==0)
				{
					bind1_1[i] = -1;
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
								bind3[g*Bp1[0]*maxl_r*maxl_r+i*maxl_r*maxl_r+ind_l2[0]*maxl_r+ind_l3[g]] = obs_r[j];
								ind_l3[g]++;
								flag=1;
								tree3[i*3]=obs_r[j];
								tree3[i*3+1]=bind2[i*maxl_r+p];
								tree3[i*3+2]=bind1_1[i];
								for(q=0; q<ind_l2[1]; q++)
								{
									tree1[i*6+q*2] = bind2[1*Bp1[0]*maxl_r+i*maxl_r+q];
									tree1[i*6+q*2+1] = bind1_1[i];
								}
								for(q=0; q<ind_l2[2]; q++)
								{
									tree1[i*2] = bind2[2*Bp1[0]*maxl_r+i*maxl_r+q];
									tree1[i*2+1] = bind1_1[i];
								}
							}
						}
						if(flag==0)
						{
							bind2[i*maxl_r+p] = -1;
							for(it=0; it<g; it++)
								ind_l3[it] = 0;
							break;
						}
						
					}
				}
				if(flag==1)
				{
					for(it=0; it<num_pre_l3; it++)
						ind_l3[it] = 0;
				}
			}
			if(flag==0)
			{
				bind1_1[i] = -1;
			}
		}
		
		for(it=0; it<num_pre_l2; it++)
			ind_l2[it] = 0;
		
	}
	
	
	
    /* local tree-pruning */
  
	int *bind2_1=(int *) malloc(num_row * sizeof(int));
	p=0;
	for(i=0; i<Bp1[0]; i++)
	{
		if( (tree3[i*3+1]!=-1) )
		{
			
			bind2_1[tree3[i*3+1]]=1;	
		}
	}
	for(i=0; i<num_row; i++)
	{
		if( (bind2_1[i]==1)&&(bind1_2[i]==1) )
			p++;
	}


	return 0;
}

