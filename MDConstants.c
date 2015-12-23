#include <stdbool.h> //bool
#include <stdio.h> //fprintf(); fscan(); sprintf();



#include "MDConstants.h"

int
Initialize (MDConstants* K, char argv[])
 {
    if (!argv)
    { 
        fprintf(stderr, "\nWrong string passed: %s;  Program aborted!!!!\n\n", argv);
        return -1;

    }

	char fname[60];

//READ INPUT.DAT
	sprintf(work_dir, "%s", argv);

	sprintf(fname, "%s/input.dat", work_dir);
	FILE *fp;
	if (!(fp = fopen(fname, "r")) )
	{
		fprintf(stderr, "\nError opening file %s;  Program aborted!!!!\n\n", fname);
		return -1;
	}
	
	
	unsigned count_scan = 0;
//TODO compress into 1 for loop
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


	if(count_scan != 17)
	{
		fprintf (stderr, "\nWrong number of input parameters, please check %s.\n", fname);
		return -1;
	}
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

	if(invalid_input)
	{
		fprintf(stderr, "\nInput parameters are not valid, please check %s.\n", fname);
		return -1;
	}

	fclose(fp);
    fp = NULL; //TODO review
 
//READ INPUT.DAT ENDS

//INIT CONSTS
    K->PartNum = 1;
    for (unsigned j = 0; j < kDIM; ++j)
    {
        K->PartNum *= K->N[j];
        K->side_minus1[j] = 1.0 / K->L[j]; 
    }
//TODO better exceptions: at the moment free runs also on non yet allocated arrays => causes segfault
    K->PartNum = (const unsigned) K->PartNum;
	K->SnapshotNum = (const unsigned) ((K->iteration_num / K->t_gap) + 1);//we count snapshot 0 as well

//INIT CONSTS ends
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
