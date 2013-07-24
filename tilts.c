#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "RC.h"
#include "flp.h"
#include "tilts.h"
extern double A[4096],B[4096];
double temp2[1024],temp3[1024];
extern double tmpr[40+EXTRA],pwr_vec[40+EXTRA]; 
extern flp_t *flp;
void generate_AB(int nodes,int inputs,double sample_intvl){
  int i,j;
  for(i=0;i<inputs;i++){
      for(j=0;j<nodes;j++){
	tmpr[j]=T_INIT+273.15;
     	pwr_vec[j]=0;
	//pwr_vec[i+20]=0;
      }
      pwr_vec[i]=1;
        
      compute_temp(pwr_vec,tmpr,flp->n_units,sample_intvl);
      for(j=0;j<nodes;j++){
	B[j*inputs+i]=tmpr[j]-T_INIT-273.15;
      }
  }

    //generate A:
    for(i=0;i<nodes;i++){
      for(j=0;j<nodes;j++){
	tmpr[j]=T_INIT+273.15;
      	pwr_vec[j]=0;
      }
      tmpr[i]+=1;
        
    //send power value into hotspot power vector
      compute_temp(pwr_vec,tmpr,flp->n_units,sample_intvl);
    
    //tilts test
      for(j=0;j<nodes;j++){
	A[j*nodes+i]=tmpr[j]-T_INIT-273.15;
      }
    }
}


void matrix_mul(double * a, double *b, double *c, int m, int n, int l)
{
	int i, j, k;
	double sum;
	for(i=0; i<m; i++) {
		for(j=0; j<l; j++) {
			sum = 0;
			for(k=0; k<n; k++) {
				sum += a[i*n+k] * b[k*l+j];
				
				
			}
			//printf("\n");
			
			c[i*l+j] = sum;
		}
	}
}

void compute_temp_matrix(double *power, double *temp, double *A, double * B, int n_input, int n_node)
{
	int i;
	//double temp2[N_NODE]; double temp3[N_NODE];
	matrix_mul((double*)A, temp, temp2, n_node, n_node, 1); 
	matrix_mul((double*)B, power, temp3, n_node, n_input, 1);
	for(i=0; i<n_node; i++) {
		temp[i] = temp2[i] + temp3[i];
	}
}

void compute_temp_tilts(double *tmpr,double *pwr,double sample_intvl,int nodes,int inputs){
  int i;
  for(i=0;i<nodes;i++){
    tmpr[i]-=(T_INIT+273.15);
  }
  compute_temp_matrix(pwr,tmpr,A,B,inputs,nodes);
   for(i=0;i<nodes;i++){
    tmpr[i]+=(T_INIT+273.15);
  }
}

void calculate_2timeAB(double *A, double *B, int n_node, int n_input)
{
	int i;
	double *A2, *B2;
	A2 = (double*)malloc(n_node * n_node * sizeof(double));
	B2 = (double*)malloc(n_node * n_input * sizeof(double));
	if(!A2 || !B2) {	//no enough memory!
		return;
	}	
	
	// Calculate A2 = A * A
	matrix_mul(A, A, A2, n_node, n_node, n_node);
	
	// calculate A + I
	for(i=0; i<n_node; i++) {
		A[i*n_node+i] += 1.0;
	}
	
	// calculate B2 = ( A + I ) * B
	matrix_mul(A, B, B2, n_node, n_node, n_input);
	
	// Copy the new matrices back into A and B.
	memcpy(A, A2, n_node * n_node * sizeof(double));
	memcpy(B, B2, n_node * n_input * sizeof(double));
	
	if(A2) free(A2);
	if(B2) free(B2);
}
