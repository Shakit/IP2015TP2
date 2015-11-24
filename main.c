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
#include <math.h>
#include <assert.h>

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
	int check;
	int solveFlag;
	glp_prob *prob;
	struct node* father;
	struct node* leftSon;
	struct node* rightSon;
} node;


/* Display node function */
void displayNode(node* n)
{
	int i;

	if (n->solveFlag == 0 && n->z == 0)
	{
		printf("NOT FEASIBLE\n");
	}
	else
	{
		printf(" z = %f\n",n->z);
		printf("Solution :\n");
		for (i = 0; i < glp_get_num_cols(n->prob); ++i) printf(" x%d = %f  \n", i+1, n->x[i]);
	}
}

/* Node creating function */
void create_node(node* n, glp_prob* prob, node* father, int y, double valy)
{
	n->father = father;
	n->leftSon = NULL;
	n->rightSon = NULL;
	n->check = 0;
	
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

	glp_write_lp(prob, NULL, "ULS.lp");

	n->solveFlag = glp_simplex(n->prob, &parm); glp_intopt(n->prob, &parmip);

	n->z = glp_mip_obj_val(n->prob);
	n->x = (double *) malloc (glp_get_num_cols(n->prob) * sizeof(double));
	for (i = 0; i < glp_get_num_cols(n->prob); ++i) n->x[i] = glp_mip_col_val(n->prob, i+1);

	
	printf("solveFlag = %d\n",n->solveFlag);
	
//	displayNode(n);
}



/* Node Checking function */
node* checked(node* n)
{
	if (n == NULL)
	{
		return NULL;
	}
	if( n->check ==0 )
	{
		return n;
	}
	else
	{	
		if( n->leftSon != NULL && n->rightSon != NULL)
		{
			node* leftSon = checked(n->leftSon);
			node* rightSon= checked(n->rightSon);

			if( leftSon != NULL || rightSon != NULL)
			{
				if (rightSon == NULL)  
				{
					return leftSon;
				}
				else {
					if (leftSon == NULL)
					{
						return rightSon;
					}
					else
					{
						if (leftSon->z <= rightSon->z)
						{
							return leftSon;
						}
						else
						{
							return rightSon;
						}
						
					}
				}
			}
			else {return NULL;}
		}
		else
		{
			return NULL;
		}
	}
}

/* base solution construction function */
node* construction (glp_prob * prob)
{
	node* res = (node *) malloc (sizeof(node));
	
	glp_prob * constProb = glp_create_prob(); glp_copy_prob(constProb, prob, GLP_ON);
	int i = glp_add_rows(constProb, 1);
	int k = glp_add_rows(constProb, 1);
    int nbj = glp_get_num_cols(prob)/3 -1;

	int ind[nbj+2];
	double val[nbj+2];
	int indk[nbj+2];
	double valk[nbj+2];
	int j;
	
	ind[0] = 0; val[0] = 1;
	ind[1] = 1; val[1] = 1;
	indk[1] = 2 ; valk[1] = 1;
	for (j = 1; j <= nbj; ++j)
	{
		ind[j] = j * 3 +2;
		indk[j] = j *3 + 3;
		val[j] = 1;
		valk[j] = 1;

		
	}

	
	glp_set_mat_row(constProb, i, nbj, ind, val);
	glp_set_row_bnds(constProb, i, GLP_FX, 0, 0);
	glp_set_mat_row(constProb, k, nbj, indk, valk);
	glp_set_row_bnds(constProb, k, GLP_FX, nbj, nbj);


	create_node(res, constProb, NULL, 0, 0);
	return res;
}

/* function that check if all Y are integer*/
int allYinteger (double* x, int size)
{
	int res = -1;
	double currentDist = 1;
	int i;
	
	for (i = 0; i < size/3; ++i)
	{
		if(x[i*3+2] != 1 && x[i*3+2] != 0)
		{
			double dist = fmin(1-x[i*3+2], x[i*3+2]);
			if(dist < currentDist)
			{
				currentDist = dist;
				res = i*3+3;
			}
		}
	}

	return res; 
}

/* Branch and bound algorithm */
node* branchAndBound (glp_prob * prob)
{
	node* root = (node *) malloc (sizeof(node));
	create_node(root, prob, NULL, 0, 0);

	node* res = construction(prob);
	double vUp = res->z; 
	
	node* node_ptr = checked(root);
	printf("------%f\n", vUp);
	displayNode(node_ptr);
	printf("------%f\n", vUp);
	int cpt = 1;
	while (node_ptr != NULL)
	{
		printf("%d\n", cpt);
		if (node_ptr->z == 0)
		{
			printf("NOT FEASIBLE!\n");
		}
		else {
			if (!(node_ptr->z >= vUp))
			{
				int y = allYinteger(node_ptr->x, glp_get_num_cols(node_ptr->prob));
				if ( y == -1)
				{
					vUp = node_ptr->z;
					res = node_ptr;
				}
				else
				{
					node* left = (node *) malloc (sizeof(node));
					node* right = (node *) malloc (sizeof(node));

					create_node(left, prob, node_ptr, y, 0);
					create_node(right, prob, node_ptr, y, 1);

					
					node_ptr->leftSon = left;
					node_ptr->rightSon = right;
				}
			}
		}

		cpt++;
	
		node_ptr->check = 1;
		node_ptr = checked(root);
	    
		printf("------------------\n");
		//	displayNode(node_ptr);
		printf("----------------------\n");
	}

	return res;
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
//	glp_simplex(prob, &parm); glp_intopt(prob, &parmip);

	/* Branch and bound */
	node* res = branchAndBound(prob);
	printf("fin BB");
	displayNode(res);

	/* Display */
//	z = glp_mip_obj_val(prob);
//	x = (double *) malloc (3*(d.nbjour+1) * sizeof(double));
//	for (i = 0; i < 3*(d.nbjour+1); ++i) x[i] = glp_mip_col_val(prob, i+1);
}
