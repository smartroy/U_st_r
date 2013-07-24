#ifndef __RC_H_
#define __RC_H_

#include "flp.h"

/* number of extra nodes due to the model	*/
/* 1 chip bottom, 5 spreader and 5 heat sink nodes (north, south, east, west and bottom)	*/
#define EXTRA	11
#define	CHIP_B	0
#define	SP_W	1
#define	SP_E	2
#define	SP_N	3
#define	SP_S	4
#define	SP_B	5
#define	SINK_W	6
#define	SINK_E	7
#define	SINK_N	8
#define	SINK_S	9
#define	SINK_B	10

/* physical constants	*/
#define RHO_SI	0.01	/* thermal resistivity of silicon between 300K-400K in (mK)/W	*/
#define	RHO_CU	0.0025	/* thermal resistivity of copper between 300K-400K in (mK)/W	*/
#define K_SI	(1.0/RHO_SI)	/* thermal conductivity of silicon	*/
#define K_CU	(1.0/RHO_CU)	/* thermal conductivity of copper	*/
#define SPEC_HEAT_SI	1.75e6	/* specfic heat of silicon in J/(m^3K)	*/
#define SPEC_HEAT_CU	3.55e6	/* specific heat of copper in J/(m^3K)	*/

/* model specific constants	*/
#define C_FACTOR	0.4		/* fitting factor to match floworks (due to lumping)	*/

#define T_INIT 40
/* set the following globals as needed	*/

/* chip specs	*/
extern double t_chip;	/* chip thickness in meters	*/
extern double thermal_threshold;	/* temperature threshold for DTM (Kelvin)*/

/* heat sink specs	*/
extern double c_convec;	/* convection capacitance in J/K */
extern double r_convec;	/* convection resistance in K/W	*/
extern double s_sink;	/* heatsink side in meters	*/
extern double t_sink;	/* heatsink thickness in meters	*/

/* heat spreader specs	*/
extern double s_spreader;	/* spreader side in meters	*/
extern double t_spreader;	/* spreader thickness in meters	*/

/* ambient temperature in kelvin	*/
extern double ambient;

/* hotspot main interface - temperature.c	*/
extern void create_RC_matrices(flp_t *flp, int omit_lateral);
extern void steady_state_temp(double *power, double *temp, int n_units);
extern void compute_temp(double *power, double *temp, int n_units, double time_elapsed);
/* differs from 'vector()' in that memory for internal nodes is also allocated	*/
extern double *hotspot_vector(int n_units);
extern void cleanup_hotspot(int n_units);
extern void set_temp (double *temp, int n_units, double val);
extern void dump_temp (flp_t *flp, double *temp, char *file);
extern void read_temp (flp_t *flp, double *temp, char *file, int clip);


/* other functions used by the above interfaces	*/

/* LUP decomposition	*/
void lupdcmp(double**a, int n, int *p);

/* get the thermal resistance values	*/
double getr(double k, double Wb, double Lb, double Ws);

/* LU forward and backward substitution	*/
void lusolve(double **a, int n, int *p, double *b, double *x);

/* 4th order Runge Kutta solver	*/
void rk4(double **c, double *y, double *pow, int n, double h, double *yout);

/* vector routines	*/
double 	*vector(int n);
void free_vector(double *v);
void dump_vector(double *v, int n);

int *ivector(int n);
void free_ivector(int *v);
void dump_ivector(int *v, int n);

/* matrix routines	*/
double **matrix(int nr, int nc);
void free_matrix(double **m, int nr);
void dump_matrix(double **m, int nr, int nc);
void copy_matrix(double **dst, double **src, int nr, int nc);

int **imatrix(int nr, int nc);
void free_imatrix(int **m, int nr);
void dump_imatrix(int **m, int nr, int nc);

void matmult(double **c, double **a, double **b, int n);
void matvectmult(double *vout, double **m, double *vin, int n);
void matinv(double **INV, double **M, int n);

#endif
