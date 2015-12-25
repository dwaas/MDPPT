#include <assert.h> //assert
#include <stdio.h> // FILE; fprintf(), fread()

#include "MDConstants.h" //MDConstants
#include "MDLoad.h" 
#include "Molecule.h" // Molecule
#include "Turbulence.h" // TurbConsts; TurbConstVecs;



//TODO const pointers vs copies
int
TurbConstsLoad 
(
    const MDConstants K,
 	TurbConsts* turb,
 	TurbConstVecs* turb_vecs
)
{
	char fname[60];	
    sprintf(fname, "%s/turbulence.pos", work_dir);
   	FILE* fp = fopen(fname, "r");

    unsigned count_scan = 0;

    for(unsigned i = 0; i < K.NF; ++i)
	{
		count_scan += fread(&turb[i].k, sizeof(double), 1, fp);
		count_scan += fread(&turb[i].omega,  sizeof(double), 1, fp);
		count_scan += fread(&turb[i].delk, sizeof(double), 1, fp);
		count_scan += fread(&turb[i].cabs, sizeof(double), 1, fp);
		count_scan += fread (&turb_vecs[i].kn, sizeof(double), kDIM, fp);
		count_scan += fread (&turb_vecs[i].c1n, sizeof(double), kDIM, fp);
		count_scan += fread (&turb_vecs[i].c2n, sizeof(double), kDIM, fp);
//TODO fread vs fscanf
	}
    fclose(fp);

    if(count_scan != K.NF * 13) goto wrong_pos_number; //TOOD macro to get the third argument of fread

   	fprintf (stderr, "\n%s read.\n", fname);
	return 0;

	wrong_pos_number:
	{
        fprintf (stderr, "Wrong number of turbulence points in %s.\n", fname);
		return -1;
    }
}



int
MDLoad 
(
    const MDConstants K, 
    Molecule* molecule, 
	TurbField* turb_field,
    const unsigned next_loop //current timestep
)
{
	//TODO check that K is initialised; check input is valid.
	assert (next_loop <= K.iteration_num);


    char fname[60];
    FILE *fp;

	sprintf(fname,"%s/t-%d.pos",work_dir, next_loop);
 	if ( !(fp = fopen(fname, "r")) ) goto wrong_file_name;
		
    unsigned count_scan = 0;
//TODO find tests for directions and turb velocities too

    double pos[kDIM];//temp variable
    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        //TODO disk I/O checks

        for (unsigned j = 0; j < kDIM; ++j)
        {
            count_scan += fread(&pos[j], sizeof(double), 1, fp);
            if ( pos[j] < K.L[j] / (-2.0) || pos[j] > K.L[j] / 2.0 )
            {
                printf("\nThe particles in the %u-th positions \
                        are not in the expexpected range for this simulation,\
                        please check %s.\n", j, fname);
                goto error;
            }//pos[] is loaded with values within the valid range

            molecule[i].position[j] = K.side_minus1[j] * pos[j];

            count_scan += fread(&molecule[i].direction[j], sizeof(double), 1, fp);
        }// molecule.position[] and molecule.direction[] are initialized		

    } 


    fclose(fp);

	if(count_scan != (K.PartNum * kDIM * 2) ) goto wrong_num_pos_entries;

	fprintf (stderr, "\n%s read.\n", fname);
	
	sprintf(fname,"%s/turbField-%d.pos",work_dir, next_loop);
    fp = fopen(fname, "r");
 	if ( !(fp = fopen(fname, "r") ) ) goto wrong_file_name;
	
    count_scan = 0;

    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        count_scan += fread(&turb_field[i], sizeof(double), kDIM, fp);
    } 

    fclose(fp);
	if(count_scan != (K.PartNum * kDIM) ) goto wrong_parameter_number;

	fprintf (stderr, "\n%s read.\n", fname);
	return 0;
//EXCEPTIONS
    wrong_parameter_number:
    {
    	fprintf (stderr, "\nFile does not contain 3 parameters, please check %s.\n", fname);
		fprintf (stderr, "\nParameters contained =  %u \n", count_scan);
        goto error; 
    }
    
    wrong_file_name:
    {
        fprintf(stderr, "\nError opening file %s\n", fname);
        goto error;
	}
    
    wrong_num_pos_entries:
	{
		printf("\nFile does not contain 6 parameters (3 positions and 3 directions), please check %s.\n", fname);
		printf("\nParameters contained =  %u \n", count_scan);
        goto error;
	}

    error:
    {
        fprintf(stderr, "\nProgram aborted!!!!\n\n");
        return -1;
    }
}

