#include "RC.h"
#include "flp.h"
#include "util.h"

#include <stdio.h>
#include <math.h>
#include <strings.h>

/* chip specs   */
double t_chip = 0.0005; /* chip thickness in meters     */
double thermal_threshold = 100 + 273.15;       /* temperature threshold for DTM (Kelvin)*/
                                                                                                                                          
/* heat sink specs      */
double c_convec = 140.4;/* convection capacitance - 140.4 J/K */
double r_convec = 0.8;  /* convection resistance - 0.8 K/W      */
double s_sink = 0.06;   /* heatsink side - 60 mm        */
double t_sink = 0.0069; /* heatsink thickness  - 6.9 mm */
                                                                                                                                          
                                                                                                                                          
/* heat spreader specs  */
double s_spreader = 0.03;       /* spreader side - 30 mm        */
double t_spreader = 0.001;      /* spreader thickness - 1 mm */

/* ambient temp in kelvin       */
double ambient = T_INIT + 273.15;   /* 35 C in kelvin       */

/* model specific globals	*/
static double factor_pack = C_FACTOR;	/* thermal capacitance fitting factor for package 	*/
static double factor_chip;				/* thermal capacitance fitting factor for silicon	*/

static double **b, **c, **inva, **invb, max_slope;

/* creates 3 matrices: invA, B, C: dT + A^-1*BT = A^-1*Power, 
 * C = A^-1 * B. note that A is a diagonal matrix (no lateral
 * capacitances. all capacitances are to ground). so, inva[i][i]
 * (= 1/a[i][i]) is just enough.
 *
 * NOTE: EXTRA nodes: 1 chip bottom, 5 spreader and 5 heat sink nodes
 * (north, south, east, west and bottom).
 */
void create_RC_matrices(flp_t *flp, int omit_lateral)
{
	int i, j, k = 0, n = flp->n_units;
	int **border;
	double **len, *gx, *gy, **g, *c_ver, **t;
	double r_sp1, r_sp2, r_hs;	/* lateral resistances to spreader and heatsink	*/

	/* NOTE: *_mid - the vertical R/C from center nodes of spreader 
	 * and heatsink. *_ver - the vertical R/C from peripheral (n,s,e,w) nodes
	 */
	double r_sp_mid, r_sp_ver, r_hs_mid, r_hs_ver, c_sp_mid, c_sp_ver, c_hs_mid, c_hs_ver;
	double gn=0, gs=0, ge=0, gw=0;
	double w_chip = get_total_width (flp);	/* x-axis	*/
	double l_chip = get_total_height (flp);	/* y-axis	*/
	FILE *fp_b,*fp_c,*fp_inva,*fp_invb;
	fp_b=fopen("B","w");
	fp_c=fopen("C","w");
	fp_invb=fopen("invB","w");
	fp_inva=fopen("invA","w");

	border = imatrix(n, 4);
	len = matrix(n, n);		/* len[i][j] = length of shared edge bet. i & j	*/
	gx = vector(n);			/* lumped conductances in x direction	*/
	gy = vector(n);			/* lumped conductances in y direction	*/
	g = matrix(n+EXTRA, n+EXTRA);	/* g[i][j] = conductance bet. nodes i & j */
	c_ver = vector(n+EXTRA);	/* vertical capacitance	*/

	b = matrix(n+EXTRA, n+EXTRA);	/* B, C, INVA  and INVB are (n+EXTRA)x(n+EXTRA) matrices	*/
	c = matrix(n+EXTRA, n+EXTRA);
	inva = matrix(n+EXTRA, n+EXTRA);
	invb = matrix(n+EXTRA, n+EXTRA);
	t = matrix (n+EXTRA, n+EXTRA);	/* copy of B	*/

	/* compute the silicon fitting factor - see pg 10 of the UVA CS tech report - CS-TR-2003-08	*/
	factor_chip = C_FACTOR * ((SPEC_HEAT_CU / SPEC_HEAT_SI) * (w_chip + 0.88 * t_spreader) \
				* (l_chip + 0.88 * t_spreader) * t_spreader / ( w_chip * l_chip * t_chip) + 1);

	/* gx's and gy's of blocks	*/
	for (i = 0; i < n; i++) {
		gx[i] = 1.0/getr(K_SI, flp->units[i].height, flp->units[i].width, l_chip);
		gy[i] = 1.0/getr(K_SI, flp->units[i].width, flp->units[i].height, w_chip);
	}

	/* shared lengths between blocks	*/
	for (i = 0; i < n; i++) 
		for (j = i; j < n; j++) 
			len[i][j] = len[j][i] = get_shared_len(flp, i, j);

	/* lateral R's of spreader and sink */
	r_sp1 = getr(K_CU, (s_spreader+3*w_chip)/4.0, (s_spreader-w_chip)/4.0, w_chip);
	r_sp2 = getr(K_CU, (3*s_spreader+w_chip)/4.0, (s_spreader-w_chip)/4.0, (s_spreader+3*w_chip)/4.0);
	r_hs = getr(K_CU, (s_sink+3*s_spreader)/4.0, (s_sink-s_spreader)/4.0, s_spreader);

	/* vertical R's and C's of spreader and sink */
	r_sp_mid = RHO_CU * t_spreader / (w_chip * l_chip);
	c_sp_mid = factor_pack * SPEC_HEAT_CU * t_spreader * (w_chip * l_chip);
	r_sp_ver = RHO_CU * t_spreader * 4.0 / (s_spreader * s_spreader - w_chip*l_chip);
	c_sp_ver = factor_pack * SPEC_HEAT_CU * t_spreader * (s_spreader * s_spreader - w_chip*l_chip) / 4.0;
	r_hs_mid = RHO_CU * t_sink / (s_spreader*s_spreader);
	c_hs_mid = factor_pack * SPEC_HEAT_CU * t_sink * (s_spreader * s_spreader);
	r_hs_ver = RHO_CU * t_sink * 4.0 / (s_sink * s_sink - s_spreader*s_spreader);
	c_hs_ver = factor_pack * SPEC_HEAT_CU * t_sink * (s_sink * s_sink - s_spreader*s_spreader) / 4.0;

	/* short the R's from block centers to a particular chip edge	*/
	for (i = 0; i < n; i++) {
		if (eq(flp->units[i].bottomy + flp->units[i].height, l_chip)) {
			gn += gy[i];
			border[i][2] = 1;	/* block is on northern border 	*/
		}	
		if (eq(flp->units[i].bottomy, 0)) {
			gs += gy[i];
			border[i][3] = 1;	/* block is on southern border	*/
		}	
		if (eq(flp->units[i].leftx + flp->units[i].width, w_chip)) {
			ge += gx[i];
			border[i][1] = 1;	/* block is on eastern border	*/
		}	
		if (eq(flp->units[i].leftx, 0)) {
			gw += gx[i];
			border[i][0] = 1;	/* block is on western border	*/
		}	
	}

	/* overall R and C between nodes */
	for (i = 0; i < n; i++) {

		/* amongst functional units	*/
		for (j = 0; j < n; j++) {
			double part = 0;
			if (!omit_lateral) {
				if (is_horiz_adj(flp, i, j)){ 
					part = gx[i] / flp->units[i].height;
					printf("%d %d horiz adj\n",i,j);
				}
				else if (is_vert_adj(flp, i,j)) {
					part = gy[i] / flp->units[i].width;
					printf("%d %d vert adj\n",i,j);
				}
			}
			g[i][j] = part * len[i][j];
		}

		/* C's from functional units to ground	*/
		c_ver[i] = factor_chip * SPEC_HEAT_SI * t_chip * flp->units[i].height * flp->units[i].width;

		/* lateral g's from block center to peripheral (n,s,e,w) spreader nodes	*/
		g[i][n+SP_N]=g[n+SP_N][i]=2.0*border[i][2]/((1.0/gy[i])+r_sp1*gn/gy[i]);
		g[i][n+SP_S]=g[n+SP_S][i]=2.0*border[i][3]/((1.0/gy[i])+r_sp1*gs/gy[i]);
		g[i][n+SP_E]=g[n+SP_E][i]=2.0*border[i][1]/((1.0/gx[i])+r_sp1*ge/gx[i]);
		g[i][n+SP_W]=g[n+SP_W][i]=2.0*border[i][0]/((1.0/gx[i])+r_sp1*gw/gx[i]);

 		/* vertical g's from block center to chip bottom */
		g[i][n+CHIP_B]=g[n+CHIP_B][i]=2.0/(RHO_SI * t_chip / (flp->units[i].height * flp->units[i].width));

	}

	/* max slope (1/vertical RC time constant) for silicon	*/
	max_slope = 1.0 / (factor_chip * t_chip * t_chip * RHO_SI * SPEC_HEAT_SI);

	/* vertical g's and C's between central nodes	*/
 	/* between chip bottom and spreader bottom */
	g[n+CHIP_B][n+SP_B]=g[n+SP_B][n+CHIP_B]=2.0/r_sp_mid;
 	/* from chip bottom to ground	*/
	c_ver[n+CHIP_B]=c_sp_mid;
 	/* between spreader bottom and sink bottom	*/
	g[n+SINK_B][n+SP_B]=g[n+SP_B][n+SINK_B]=2.0/r_hs_mid;
 	/* from spreader bottom to ground	*/
	c_ver[n+SP_B]=c_hs_mid;
 	/* from sink bottom to ground	*/
	c_ver[n+SINK_B]=c_convec;

	/* g's and C's from peripheral(n,s,e,w) nodes	*/
	for (i = 1; i <= 4; i++) {
 		/* vertical g's between peripheral spreader nodes and spreader bottom */
		g[n+SP_B-i][n+SP_B]=g[n+SP_B][n+SP_B-i]=2.0/r_sp_ver;
 		/* lateral g's between peripheral spreader nodes and peripheral sink nodes	*/
		g[n+SP_B-i][n+SINK_B-i]=g[n+SINK_B-i][n+SP_B-i]=2.0/(r_hs + r_sp2);
 		/* vertical g's between peripheral sink nodes and sink bottom	*/
		g[n+SINK_B-i][n+SINK_B]=g[n+SINK_B][n+SINK_B-i]=2.0/r_hs_ver;
 		/* from peripheral spreader nodes to ground	*/
		c_ver[n+SP_B-i]=c_sp_ver;
 		/* from peripheral sink nodes to ground	*/
		c_ver[n+SINK_B-i]=c_hs_ver;
	}

	/* calculate matrices A, B such that A(dT) + BT = POWER */

	for (i = 0; i < n+EXTRA; i++) {
		for (j = 0; j < n+EXTRA; j++) {
			if (i==j) {
				inva[i][j] = 1.0/c_ver[i];
				if (i == n+SINK_B)	/* sink bottom */
					b[i][j] += 1.0 / r_convec;
				for (k = 0; k < n+EXTRA; k++) {
					if ((g[i][k]==0.0)||(g[k][i])==0.0) 
						continue;
					else 
					/* here is why the 2.0 factor comes when calculating g[][]	*/
						b[i][j] += 1.0/((1.0/g[i][k])+(1.0/g[k][i]));
				}
			} else {
				inva[i][j]=0.0;
				if ((g[i][j]==0.0)||(g[j][i])==0.0)
					b[i][j]=0.0;
				else
				b[i][j]=-1.0/((1.0/g[i][j])+(1.0/g[j][i]));
			}
		}
	}

	/* we are always going to use the eqn dT + A^-1 * B T = A^-1 * POWER. so, store  C = A^-1 * B	*/
	matmult(c, inva, b, n+EXTRA);
	/* we will also be needing INVB so store it too	*/
	copy_matrix(t, b, n+EXTRA, n+EXTRA);
	matinv(invb, t, n+EXTRA);
	for (i = 0; i < n+EXTRA; i++) {
		for (j = 0; j < n+EXTRA; j++) {
			fprintf(fp_inva,"%f  ",inva[i][j]);
			fprintf(fp_invb,"%f  ",invb[i][j]);
			fprintf(fp_c,"%f  ",c[i][j]);
			fprintf(fp_b,"%f  ",b[i][j]);
		}
		fprintf(fp_inva, "\n");
		fprintf(fp_invb, "\n");
		fprintf(fp_c , "\n");
		fprintf(fp_b, "\n");
	}
	fclose(fp_inva);
	fclose(fp_b);
	fclose(fp_c);
	fclose(fp_invb);

/*	dump_vector(c_ver, n+EXTRA);	*/
/*	dump_matrix(invb, n+EXTRA, n+EXTRA);	*/
/*	dump_matrix(c, n+EXTRA, n+EXTRA);	*/

	/* cleanup */
	free_matrix(t, n+EXTRA);
	free_matrix(g, n+EXTRA);
	free_matrix(len, n);
	free_imatrix(border, n);
	free_vector(c_ver);
	free_vector(gx);
	free_vector(gy);
}

/* setting internal node power numbers	*/
void set_internal_power (double *pow, int n_units)
{
	int i;
	for (i=n_units; i < n_units+SINK_B; i++)
		pow[i] = 0;
	pow[n_units+SINK_B] = ambient / r_convec;
}

/* power and temp should both be alloced using hotspot_vector. 
 * 'b' is the 'thermal conductance' matrix. i.e, b * temp = power
 *  => temp = invb * power
 */
void steady_state_temp(double *power, double *temp, int n_units) 
{
	/* set power numbers for the virtual nodes */
	set_internal_power(power, n_units);

	/* find temperatures	*/
	matvectmult(temp, invb, power, n_units+EXTRA);
}

/* required precision in degrees	*/
#define PRECISION	0.1
#define TOO_LONG	100000
#define MIN_ITER	1

/* compute_temp: solve for temperature from the equation dT + CT = inv_A * Power 
 * Given the temperature (temp) at time t, the power dissipation per cycle during the 
 * last interval (time_elapsed), find the new temperature at time t+time_elapsed.
 * power and temp should both be alloced using hotspot_vector
 */

void compute_temp(double *power, double *temp, int n_units, double time_elapsed)
{
	int i;
	double *pow, h, n_iter;
	pow = vector(n_units+EXTRA);
	
	/* set power numbers for the virtual nodes */
	set_internal_power(power, n_units);

	/* find (inv_A)*POWER */
	matvectmult(pow, inva, power, n_units+EXTRA);

 	/* step size for 4th-order Runge-Kutta  - assume worst case	*/
	h = PRECISION / max_slope;
	n_iter = time_elapsed / h;
	n_iter = (n_iter > MIN_ITER) ? n_iter : MIN_ITER;	/* do atleast MIN_ITER iterations	*/
	h = time_elapsed / n_iter;
	
	if (n_iter >= TOO_LONG)
		fprintf(stderr, "warning: calling interval too large, performing %.0f iterations - it may take REALLY long\n", n_iter);
	
	/* Obtain temp at time (t+h). 
	 * Instead of getting the temperature at t+h directly, we do it 
	 * in n_iter steps to reduce the error due to rk4
	 */
	for (i = 0; i < n_iter; i++) 	
		rk4(c, temp, pow, n_units+EXTRA, h, temp);

	free_vector(pow);
}
void compute_temp_simple(double *power,double *temp,double interval){
	int i;
	double thermal_r,thermal_c;
	thermal_r=5.6;thermal_c=0.061;
	temp[0]=T_INIT+273.15+power[0]*thermal_r-(power[0]*thermal_r-temp[0]+T_INIT+273.15)*0.746;//exp(-interval/(thermal_r*thermal_c));
	temp[1]=T_INIT+273.15+power[1]*thermal_r-(power[1]*thermal_r-temp[1]+T_INIT+273.15)*0.746;//exp(-interval/(thermal_r*thermal_c));

}

/* differs from 'vector()' in that memory for internal nodes is also allocated	*/
double *hotspot_vector(int n_units)
{
	return vector(n_units+EXTRA);
}

/* sets the temperature of a vector 'temp' allocated using 'hotspot_vector'	*/
void set_temp(double *temp, int n_units, double val)
{
	int i;
	for(i=0; i < n_units + EXTRA; i++)
		temp[i] = val;
}

/* dump temperature vector alloced using 'hotspot_vector' to 'file' */ 
void dump_temp(flp_t *flp, double *temp, char *file)
{
	int i;
	char str[STR_SIZE];
	FILE *fp = fopen (file, "w");
	if (!fp) {
		sprintf (str,"error: %s could not be opened for writing\n", file);
		fatal(str);
	}
	/* on chip temperatures	*/
	for (i=0; i < flp->n_units; i++)
		fprintf(fp, "%s\t%.1f\n", flp->units[i].name, temp[i]);

	/* internal node temperatures	*/
	for (i=0; i < EXTRA; i++) {
		sprintf(str, "inode_%d", i);
		fprintf(fp, "%s\t%.1f\n", str, temp[i+flp->n_units]);
	}
	fclose(fp);	
}

/* 
 * read temperature vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_vector'. values are clipped to thermal
 * threshold based on 'clip'
 */ 
void read_temp(flp_t *flp, double *temp, char *file, int clip)
{
	int i, idx;
	double max=0, val;
	char str[STR_SIZE], name[STR_SIZE];
	FILE *fp = fopen (file, "r");
	if (!fp) {
		sprintf (str,"error: %s could not be opened for reading\n", file);
		fatal(str);
	}	

	/* find max temp on the chip	*/
	for (i=0; i < flp->n_units; i++) {
		fgets(str, STR_SIZE, fp);
		if (feof(fp))
			fatal("not enough lines in temperature file\n");
		if (sscanf(str, "%s%lf", name, &val) != 2)
			fatal("invalid temperature file format\n");
		idx = get_blk_index(flp, name);
		if (idx >= 0)
			temp[idx] = val;
		else	/* since get_blk_index calls fatal, the line below cannot be reached	*/
			fatal ("unit in temperature file not found in floorplan\n");
		if (temp[idx] > max)
			max = temp[idx];
	}

	/* internal node temperatures	*/	
	for (i=0; i < EXTRA; i++) {
		fgets(str, STR_SIZE, fp);
		if (feof(fp))
			fatal("not enough lines in temperature file\n");
		if (sscanf(str, "%s%lf", name, &val) != 2)
			fatal("invalid temperature file format\n");
		sprintf(str, "inode_%d", i);
		if (strcasecmp(str, name))
			fatal("invalid temperature file format\n");
		temp[i+flp->n_units] = val;	
	}

	fclose(fp);	

	/* clipping	*/
	if (clip && (max > thermal_threshold)) {
		/* if max has to be brought down to thermal_threshold, 
		 * (w.r.t the ambient) what is the scale down factor?
		 */
		double factor = (thermal_threshold - ambient) / (max - ambient);
	
		/* scale down all temperature differences (from ambient) by the same factor	*/
		for (i=0; i < flp->n_units + EXTRA; i++)
			temp[i] = (temp[i]-ambient)*factor + ambient;
	}
}


void cleanup_hotspot(int n_units)
{
	free_matrix(inva, n_units+EXTRA);
	free_matrix(b, n_units+EXTRA);
	free_matrix(invb, n_units+EXTRA);
	free_matrix(c, n_units+EXTRA);
}
