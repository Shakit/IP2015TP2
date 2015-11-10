/* main.c
 * Authors : DELAVERNHE Florian, LEGRU Guillaume
 * Date 10 11 2015
 *
 * Implémentation d'un algorithme de Branch & Bound
 *
 */

/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glpk.h>

typedef struct
{
	int nbjour;
	int* d;
	int* p;
	int* h;
	int* f;
} data;

void filereader(char* filename, data* d)
{
	FILE *fin;
	int i;
	int val;

	fin = fopen(filename, "r");

	fscanf(fin, "%d", &val);
	d->nbjour = val;

	/* Allocations */
	d->d = (int *) malloc ((d->nbjour +1) * sizeof(int));	
	d->p = (int *) malloc ((d->nbjour +1) * sizeof(int));	
	d->h = (int *) malloc ((d->nbjour +1) * sizeof(int));	
	d->f = (int *) malloc ((d->nbjour +1) * sizeof(int));

	
    /* Read and fill*/	
	d->d[0] = 0;
	d->p[0] = 0;
	d->h[0] = 0;
	d->f[0] = 0;
				
	for (i = 1; i <= d->nbjour; ++i)
	{
		fscanf(fin, "%d", &val);
		d->d[i] = val;
	}
	for (i = 1; i <= d->nbjour; ++i)
	{
		fscanf(fin, "%d", &val);
		d->p[i] = val;
	}
	for (i = 1; i <= d->nbjour; ++i)
	{
		fscanf(fin, "%d", &val);
		d->h[i] = val;
	}
	for (i = 1; i <= d->nbjour; ++i)
	{
		fscanf(fin, "%d", &val);
		d->f[i] = val;
	}

	fclose(fin);
}

int main(int argc, char** argv)
{
	/*==================================================*/
	/* Variables */
	data d;
	int* M;

	/* GLPK */
	int *ia, *ja, *ar;
	double z;
	double *x;

	/* Misc */
	int i;

	/* Check up */
	if(argc != 2)
	{
		printf("ERROR : no data file !");
		exit(1);		
	}
	
	/* Initialization */
	filereader(argv[1], &d);

	if(d.nbjour < 1)
	{
		printf("Obvious...");
		return 0;	
	}
	
	M = (int*) malloc ((d.nbjour +1)* sizeof(int));
	M[d.nbjour] = d.d[nbjour];
	for (i = d.nbjour-1; i >=0; --i)
	{
		M[i] = d.d[i] + M[i+1];
	}

	/* Problem creation*/
	glp_prob *prob;
	prob = glp_create_prob();
	glp_set_prob_name(prob, "ULS");
	glp_set_obj_dir(prob, GLP_MIN);

	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF;

	glp_smcp parmip;
	glp_init_iocp(&parmip);
	parmip.msg_lev = GLP_MSG_OFF;

	/* Number of constraints : 2 * nbjour +2 */
	glp_add_rows(prob, 2*d.nbjour +2);
	for (i = 1; i <= d.nbjour; ++i)
	{
		glp_set_raw_bnds(prob, i, GLP_FX, d.d[i], d.d[i]);
	}
	for (i = d.nbjour+1; i < 2*d.nbjour; ++i)
	{
		glp_set_raw_bnds(prob, i, GLP_LO, 0, 0);
	}	
	glp_set_raw_bnds(prob, 2*d.nbjour+1, GLP_FX, 0.0, 0.0);
	glp_set_raw_bnds(prob, 2*d.nbjour+2, GLP_FX, 0.0, 0.0);

	/* Number of variables : 3*(nbjour +1)*/
	glp_add_cols(prob, 3*(d.nbjour+1));
	for (i = 0; i < d.nbjour +1; ++i)
	{
		glp_set_col_bnds(prob, i*3 +1, GLP_LO, 0.0, 0.0);
		glp_set_col_bnds(prob, i*3 +2, GLP_LO, 0.0, 0.0);
		glp_set_col_bnds(prob, i*3 +3, GLP_DB, 0.0, 1.0);
	}


		


}
