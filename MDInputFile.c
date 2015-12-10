#include <stdbool.h> //bool
#include <stdio.h> //fprintf(); fscan(); sprintf();
#include <stdlib.h> 



#include "MDInputFile.h"




void
Initialize (MDConstants* K, char argv[])
 {

	char work_dir[40], fname[60];

//READ INPUT.DAT
	sprintf(work_dir, "%s", argv);

	sprintf(fname, "%s/input.dat", work_dir);
	FILE *fp;
	if ((fp = fopen(fname, "r"))==NULL)
	{
		fprintf(stderr, "\nError opening file %s;  Program aborted!!!!\n\n", fname);
		exit(1);
	}
	
	
	unsigned count_scan = 0;

	count_scan += fscanf(fp, "%d", K->iteration_num);
	count_scan += fscanf(fp, "%d %d %d", K->Nx, K->Ny, K->Nz);
	count_scan += fscanf(fp, "%lf %lf %lf", K->Lx,  K->Ly, K->Lz );
	count_scan += fscanf(fp, "%lf", K->dcut);
	count_scan += fscanf(fp, "%s", K->starting_branch);
	count_scan += fscanf(fp, "%lf", K->v_0); /* LC alignment: constant velocity*/
	count_scan += fscanf(fp, "%lf", K->gamma_tilde); /* LC alignment: strength of perturbation*/
	count_scan += fscanf(fp, "%lf", K->delta_t); /*  time step */
	count_scan += fscanf(fp, "%lf", K->gamma_rel); /*  relaxation constant */
	count_scan += fscanf(fp, "%lf %d", K->kappa, K->NF);
	count_scan += fscanf(fp, "%d %d", K->t_gap, K->deltaS);

	bool invalid_input = K->iteration_num <= 0 ||
								K->Nx <= 0 || K->Ny <= 0 || K->Nz <= 0 ||
								K->Lx <= 0 || K->Ly <= 0 || K->Lz <= 0 ||
								K->dcut <= 0 ||
								K->v_0 < 0 || 
								K->delta_t <= 0 ||
								K->gamma_rel <= 0 ||
								K->kappa < 0 || K->NF == 0 || 
								K->t_gap <= 0 || K->deltaS <= 0;
//kappa and v_0 can be zero
//TODO nonnegative unsigned
	if(count_scan<16)
	{
		printf("\nNot enough input parameters, please check %s.\n", fname);
		exit(1);
	}

	if(invalid_input)
	{
		printf("\nInput parameters are not valid, please check %s.\n", fname);
		exit(1);
	}

	fclose(fp);
    fp = NULL;
 
//READ INPUT.DAT ENDS
//TODO fix progress  bar
//INIT CONSTS
	K->DimNum  = (const unsigned) 3; //TODO more elegant way
    K->PartNum = (const unsigned) K->Nx * K->Ny * K->Nz;
	K->Length = (const unsigned) K->Lx; // cubic box
	K->InvLength = (const double) 1./K->Length; 
	K->SnapshotNum = (const unsigned) ((K->iteration_num / K->t_gap) + 1);//we count snapshot 0 as well

//INIT CONSTS ends
	fprintf(stderr, "\niteration num = %d t_gap =%i \
						\nsnapshot number = %i\
						\nparticle num = %i length = %f\n",
			K->iteration_num, K->t_gap, 
            K->SnapshotNum, 
            K->PartNum, K->Length); //TODO give meaningful info

    return;
}


