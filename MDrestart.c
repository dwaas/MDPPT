#include <assert.h> //assert
#include <stdio.h> // FILE; fprintf(), fread()
#include <stdlib.h> //exit(); 

#include "MDInputFile.h" //MDConstants
#include "Molecule.h" // Molecule
#include "Turbulence.h" // TurbConsts; TurbConstVecs;


void
MDLoad 
(
    const MDConstants K, 
    Molecule* molecule, 
    TurbConsts* turb,
    TurbConstVecs* turb_vecs,
    const unsigned next_loop //current timestep
)
{
	//TODO check that K is initialised
    unsigned molecule_size = sizeof (molecule) / 
								sizeof (Molecule*);
	assert (molecule_size == K.PartNum); 
	unsigned turb_size = sizeof (turb) /
						 sizeof (TurbConsts*);
	assert (turb_size == K.NF); 
	unsigned turb_vecs_size = sizeof (turb_vecs) / 
								sizeof (TurbConsts*);
	assert (turb_vecs_size == K.NF); 
	assert (next_loop <= K.iteration_num);

	unsigned last_t;
    double x1, y1, z1;

    char fname[60];
    FILE *fp;
    unsigned count_scan;
    double 
            side_x_minus1 = 1.0 / K.Lx,
            side_y_minus1 = 1.0 / K.Ly,
            side_z_minus1 = 1.0 / K.Lz;

	sprintf(fname,"%s/t-%d.pos",work_dir, next_loop);
    fp = fopen(fname, "r");
 	if ((fp = fopen(fname, "r"))==NULL) //TODO refactor
		{
			fprintf(stderr, "\nError opening file %s;  Program aborted!\n\n", fname);
			exit(1);
		}

   count_scan = 0;

    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        count_scan += fread(&x1, sizeof(double), 1, fp);
        count_scan += fread(&y1, sizeof(double), 1, fp);
        count_scan += fread(&z1, sizeof(double), 1, fp);

        molecule[i].x = side_x_minus1 * x1;
        molecule[i].y = side_y_minus1 * y1; 
        molecule[i].z = side_z_minus1 * z1; 

        count_scan += fread(&molecule[i].e_x, sizeof(double), 1, fp);
        count_scan += fread(&molecule[i].e_y, sizeof(double), 1, fp);
        count_scan += fread(&molecule[i].e_z, sizeof(double), 1, fp);
    } 

//only used to reload    count_scan += fread(&last_t, sizeof(int), 1, fp);

    fclose(fp);
//TODO to refactor; make clear if evaluation
	if(count_scan != 6)
	{
		printf("\nFile does not contain 6 parameters (3 positions and 3 directions), please check %s.\n", fname);
		printf("\nParameters contained =  %u \n", count_scan);

		exit(1);
	}
	
    printf("Read restart.pos.\n");

    sprintf(fname, "%s/turbulence.pos", work_dir);
    fp = fopen(fname, "r");
    count_scan = 0;
    for(unsigned i = 0; i < K.NF; ++i)
    {
        count_scan += fread(&turb[i].k, sizeof(double), 1, fp);
        count_scan += fread(&turb[i].omega,  sizeof(double), 1, fp);
        count_scan += fread(&turb[i].delk, sizeof(double), 1, fp);
        count_scan += fread(&turb[i].cabs, sizeof(double), 1, fp);

        count_scan += fread(&turb_vecs[i].kn_x, sizeof(double), 1, fp);
        count_scan += fread(&turb_vecs[i].kn_y, sizeof(double), 1, fp);
        count_scan += fread(&turb_vecs[i].kn_z, sizeof(double), 1, fp);

        count_scan += fread(&turb_vecs[i].c1n_x, sizeof(double), 1, fp);
        count_scan += fread(&turb_vecs[i].c1n_y, sizeof(double), 1, fp);
        count_scan += fread(&turb_vecs[i].c1n_z, sizeof(double), 1, fp);

        count_scan += fread(&turb_vecs[i].c2n_x, sizeof(double), 1, fp);
        count_scan += fread(&turb_vecs[i].c2n_y, sizeof(double), 1, fp);
        count_scan += fread(&turb_vecs[i].c2n_z, sizeof(double), 1, fp);
    }
    fclose(fp);

    if(count_scan != K.NF * 13) //TODO refactor
    {
        printf("Wrong number of turbulence points in %s.\n", fname);
        exit(0);
    }
    printf("Read turbulence.pos.\n");
//TODO load TURB files
	return;	
    //return last_t+1; //loop_start
}



/* 
		
	

		}
		

		fclose(fp);
*/

