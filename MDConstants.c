#include <assert.h>
#include <stdbool.h> //bool
#include <stdio.h> //fprintf(); fscan(); sprintf();

#include "debug.h"
#include "MDConstants.h"

int
Initialize 
(
	MDConstants* K, 
	const char fname[]
)
{
	FILE *fp = NULL;
	fp = fopen(fname, "r");
	check (fp, "\nWrong filename passed: %s\n", fname);


	unsigned count_scan = 0;

	count_scan += fscanf(fp, "%u", &(K->iteration_num) );
	CountReadUnsigned (&count_scan, fp, K->N, kDIM);
	CountReadDouble (&count_scan, fp, K->L, kDIM);
	count_scan += fscanf(fp, "%lf", &(K->dcut) );
	count_scan += fscanf(fp, "%s", (K->starting_branch) );
	count_scan += fscanf(fp, "%lf", &(K->v_0) ); /* LC alignment: constant velocity*/
	count_scan += fscanf(fp, "%lf", &(K->gamma_tilde) ); /* LC alignment: strength of perturbation*/
	count_scan += fscanf(fp, "%lf", &(K->delta_t) ); /*  time step */
	count_scan += fscanf(fp, "%lf", &(K->gamma_rel) ); /*  relaxation constant */
	count_scan += fscanf(fp, "%lf %d", &(K->kappa) , &(K->NF) );
	count_scan += fscanf(fp, "%d %d", &(K->t_gap) , &(K->deltaS) );


	check ( (count_scan == 17), 
			"\nWrong number of input parameters, please check:\
			\n%s\n", 
			fname
	      );

	//TODO range evaluation 

	bool invalid_input = K->iteration_num <= 0 ||
		K->dcut <= 0 ||
		K->v_0 < 0 || 
		K->delta_t <= 0 ||
		K->gamma_rel <= 0 ||
		K->kappa < 0 || K->NF == 0 || 
		K->t_gap <= 0 || K->deltaS <= 0;

	for (unsigned j = 0; j < kDIM; ++j)
	{
		invalid_input = invalid_input || (K->N[j] <= 0);
		invalid_input = invalid_input || (K->L[j] <= 0);
	}

	//kappa and v_0 can be zero
	check (
			!invalid_input,
			"\nInvalid value ranges, please check:\
			\n%s.\n",
			fname
	      );

	fclose(fp);
	fp = NULL; //TODO review

	//READ INPUT.DAT ENDS
	CalcConsts (K); //Derived constants

	fprintf
		(
		 stderr, 
		 "\niteration num = %u t_gap = %u \
		 \nsnapshot number = %i\
		 \nparticle num = %i \n",
		 K->iteration_num, K->t_gap, 
		 K->SnapshotNum, 
		 K->PartNum
		);

	//TODO unroll dimensional loops
	PrintVals ("L_%u = %f\t", K->L, kDIM);
	PrintVals ("side_minus1_%u = %f\t", K->side_minus1, kDIM);

	return 0;

error:
	{
		fprintf(stderr, "\nProgram aborted!!!!\n\n");
		return -1;
	}
}

void
CalcConsts
(
	MDConstants* K
)
{
	//INIT CONSTS
    K->PartNum = 1;
    for (unsigned j = 0; j < kDIM; ++j)
    {
        K->PartNum *= K->N[j];
        K->side_minus1[j] = 1.0 / K->L[j]; 
    }
    K->PartNum = (const unsigned) K->PartNum;
	K->SnapshotNum = (const unsigned) ((K->iteration_num / K->t_gap) + 1);//we count snapshot 0 as well


	return;
}
void
CountReadDouble 
(
    unsigned* count_scan,
    FILE* fp,
    double* v_ptr,
    unsigned dim
)
{
    for (unsigned j = 0; j < dim; ++j, ++v_ptr)
    {
        *count_scan += fscanf(fp, "%lf", v_ptr );
    }

    return;
}
 
void
CountReadUnsigned
(
    unsigned* count_scan,
    FILE* fp,
    unsigned* v_ptr,
    const unsigned dim
)
{
    for (unsigned j = 0; j < dim; ++j, ++v_ptr)
    {
        *count_scan += fscanf(fp, "%u", v_ptr );
    }

    return;
}
 
void
PrintVals
(
    const char message[],
    double* v_ptr,
    unsigned dim
)
{
    for (unsigned j = 0; j < dim; ++j, ++v_ptr)
    {
        fprintf (stderr, message, j, *v_ptr);
    } 
    fprintf (stderr, "\n");
    return;
} 
