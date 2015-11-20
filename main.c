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

/* Data structure */
typedef struct data
{
	int nbjour;
	int* d;
	int* p;
	int* h;
	int* f;
} data;

/* Node structure */
typedef struct node
{
    double* x;
	double z;
	int checked;
	glp_prob *prob;
	struct node* father;
	struct node* leftSon;
	struct node* rightSon;
} node;

/* Node creating function */
void create_node(node* n, glp_prob* prob, node* father, int y, double valy)
{
	n->father = father;
	n->leftSon = NULL;
	n->rightSon = NULL;
	
	int i = 0;
	int ind[] = {0,y};
	double val[] = {0,1}; 

	if (n-> father == NULL)
	{
		n->prob = prob;
	}
	else
	{
		n->prob = glp_create_prob();
		glp_copy_prob(n->prob, n->father->prob, GLP_ON);
		i = glp_add_rows(n->prob, 1);
		glp_set_mat_row(n->prob, i, 1, ind, val);
		glp_set_row_bnds(n->prob, i, GLP_FX, valy, valy);
	}

	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF;

	glp_iocp parmip;
	glp_init_iocp(&parmip);
	parmip.msg_lev = GLP_MSG_OFF;

	glp_simplex(n->prob, &parm); glp_intopt(n->prob, &parmip);

	n->z = glp_mip_obj_val(n->prob);
	n->x = (double *) malloc (glp_get_num_cols(n->prob) * sizeof(double));
	for (i = 0; i < glp_get_num_cols(n->prob); ++i) n->x[i] = glp_mip_col_val(n->prob, i+1);
}

/* Display node function */
void displayNode(node* n)
{
	int i;
	
	printf(" z = %f\n",n->z);
	printf("Solution :\n");
	for (i = 0; i < glp_get_num_cols(n->prob); ++i) printf(" x%d = %f  \n", i, n->x[i]);
}

/* File reading function */
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
	int *ia, *ja;
	double *ar;
	double z;
	double *x;
	int notZeroCount;

	/* Misc */
	int i, pos;

	/* Check up */
	if(argc != 2)
	{
		printf("ERROR : no data file !\n");
		exit(1);		
	}
	
	/* Initialization */
	filereader(argv[1], &d);

	if(d.nbjour < 1)
	{
		printf("Obvious...\n");
		return 0;	
	}
	
	M = (int*) malloc ((d.nbjour +1)* sizeof(int));
	M[d.nbjour] = d.d[d.nbjour];
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

	glp_iocp parmip;
	glp_init_iocp(&parmip);
	parmip.msg_lev = GLP_MSG_OFF;

	/* Number of constraints : 2 * nbjour +2 */
	glp_add_rows(prob, 2*d.nbjour +2);
	for (i = 1; i <= d.nbjour; ++i)
	{
		glp_set_row_bnds(prob, i, GLP_FX, d.d[i], d.d[i]);
	}
	for (i = d.nbjour+1; i <= 2*d.nbjour; ++i)
	{
		glp_set_row_bnds(prob, i, GLP_LO, 0, 0);
	}	
	glp_set_row_bnds(prob, 2*d.nbjour+1, GLP_FX, 0.0, 0.0);
	glp_set_row_bnds(prob, 2*d.nbjour+2, GLP_FX, 0.0, 0.0);

	/* Number of variables : 3*(nbjour +1)*/
	glp_add_cols(prob, 3*(d.nbjour+1));
	for (i = 0; i < d.nbjour +1; ++i)
	{
		glp_set_col_bnds(prob, i*3 +1, GLP_LO, 0.0, 0.0);
		glp_set_col_bnds(prob, i*3 +2, GLP_LO, 0.0, 0.0);
		glp_set_col_bnds(prob, i*3 +3, GLP_DB, 0.0, 1.0);
	}
	for (i = 1; i <= 3*(d.nbjour+1); ++i)
	{
		glp_set_col_kind(prob, i, GLP_CV);
	}

	/* Coefficients of the economic function */
	glp_set_obj_coef(prob, 1, 0);
	glp_set_obj_coef(prob, 2, 0);
	glp_set_obj_coef(prob, 3, 0);
	for (i = 1; i <=d.nbjour; ++i)
	{
		glp_set_obj_coef(prob, 3*i+1, d.p[i]);
		glp_set_obj_coef(prob, 3*i+2, d.h[i]);
		glp_set_obj_coef(prob, 3*i+3, d.f[i]);
	}

	/* Matrix */
	notZeroCount = 5 * d.nbjour + 2;

	ia = (int *) malloc ((1+notZeroCount) * sizeof(int));
	ja = (int *) malloc ((1+notZeroCount) * sizeof(int));	
	ar = (double *) malloc ((1+notZeroCount) * sizeof(double));

	pos = 1;
	for (i = 1; i <= d.nbjour; ++i)
	{
		ia[pos] = i;
		ia[pos+1] = i;
		ia[pos+2] = i;

		ja[pos] = i*3-1;
		ja[pos+1] = i*3+1;
		ja[pos+2] = i*3+2;

		ar[pos] = 1.0;
		ar[pos+1] = 1.0;
		ar[pos+2] = -1.0;
			
		pos += 3;
	}

	for (i = 1; i <= d.nbjour; ++i)
	{
		ia[pos] = i + d.nbjour;
		ia[pos+1] = i + d.nbjour;

		ja[pos] = i*3+1;
		ja[pos+1] = i*3+3;

		ar[pos] = -1.0;
		ar[pos+1] = M[i];
			
		pos += 2;
	}
	
	ia[pos] = 2*d.nbjour +1;
	ia[pos+1] = 2 * d.nbjour+2;

	ja[pos] = 3*(d.nbjour+1)-1;
	ja[pos+1] = 2;

	ar[pos] = 1.0;
	ar[pos+1] = 1.0;

	pos += 2;

	glp_load_matrix(prob, notZeroCount, ia, ja , ar);

	/* Writing in a file */
	glp_write_lp(prob, NULL, "ULS.lp");

	/* Solve */
	glp_simplex(prob, &parm); glp_intopt(prob, &parmip);

	/*Branch and bound*/
	node* n = (node*) malloc (sizeof(node));
	create_node(n, prob, NULL, 0, 0); // first node 0
	displayNode(n);
	node* zbra = (node*) malloc (sizeof(node));
	
	create_node(zbra,prob,n,6,1);
	displayNode(zbra);
	
	/* Display */
//	z = glp_mip_obj_val(prob);
//	x = (double *) malloc (3*(d.nbjour+1) * sizeof(double));
//	for (i = 0; i < 3*(d.nbjour+1); ++i) x[i] = glp_mip_col_val(prob, i+1);
}
