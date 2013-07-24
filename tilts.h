#ifndef __TILTS_H_
#define __TILTS_H_
#include <stdlib.h>
#include <stdio.h>
#include "RC.h"
#include "flp.h"

void generate_AB(int nodes,int inputs,double sample_intvl);
void matrix_mul(double * a, double *b, double *c, int m, int n, int l);
void compute_temp_matrix(double *power, double *temp, double *A, double * B, int n_input, int n_node);
void compute_temp_tilts(double *tmpr,double *pwr,double sample_intvl,int nodes,int inputs);
void calculate_2timeAB(double *A, double *B, int n_node, int n_input);
#endif
