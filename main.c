#include <assert.h> //assert()
#include <stdbool.h> //true, false
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "debug.h"

#include "MDConstants.h" //MDConstants;
#include "MDLoad.h" // MDLoad();
#include "MDPostProcessing.h" // MeanKineticEnergy()
#include "Molecule.h" // Molecule
#include "Turbulence.h" //TurbConsts; TurbConstVecs

#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))

//TODO documentation
//TODO testsuite: wrong args, wrong input file, wrong data in input file, no memory


int
main (int argc, char *argv[])
{
	check (
			argc == 2,
			"\nNumber of arguments is wrong! \
			\nusage: ./post_proc working_directory_path!!!!\n\n"
		  );

    char work_dir[40];
	sprintf (work_dir, "%s", argv[1]);


    Molecule** positions = NULL;
    TurbField** turb_velocities = NULL;
    Tensor2** strain_rate = NULL;    

    TurbConsts* turb = NULL;
    TurbConstVecs* turb_vecs = NULL;
    KraichnanMode*** kraich_modes = NULL;


    //CONSTANTS
    MDConstants K;
    char input_dat[60];
	sprintf (input_dat, "%s/input.dat", work_dir);
	check 
		(
			!Initialize (&K, input_dat), 
			"Constant init failed" 
		);

    //allocate 
    positions = (Molecule**) calloc (K.SnapshotNum, sizeof (Molecule*) );
	check_mem (positions);

    turb_velocities = (TurbField**) calloc (K.SnapshotNum, sizeof (TurbField*) );
    if (!turb_velocities) goto no_memory;

    strain_rate = (Tensor2**) calloc (K.SnapshotNum, sizeof (Tensor2*) );
    if (!strain_rate) goto no_memory;


    kraich_modes = (KraichnanMode***) calloc (K.SnapshotNum, sizeof (KraichnanMode**) );
    if (!kraich_modes) goto no_memory;

    for (unsigned n = 0; n < K.SnapshotNum; ++n) 
    {
        positions[n] = (Molecule*) calloc (K.PartNum, sizeof (Molecule) );
        if (!positions[n]) goto no_memory;

        turb_velocities[n] = (TurbField*) calloc (K.PartNum, sizeof (TurbField) );
        if (!turb_velocities[n]) goto no_memory;

        strain_rate[n] = (Tensor2*) calloc (K.PartNum, sizeof (Tensor2) );
        if (!strain_rate[n]) goto no_memory;
		//kraich_modes[n] = NULL;
        kraich_modes[n] = (KraichnanMode**) calloc (K.PartNum, sizeof (KraichnanMode*) );
        if (!kraich_modes[n]) goto no_memory;

        for (unsigned i = 0; i < K.PartNum; ++i)
        {
            kraich_modes[n][i] = (KraichnanMode*) calloc (K.NF, sizeof (KraichnanMode) );
            if (!kraich_modes[n][i]) goto no_memory;

        }

    }

    turb = (TurbConsts*) calloc (K.NF, sizeof (TurbConsts) );
    if (!turb) goto no_memory;

    turb_vecs = (TurbConstVecs*) calloc (K.NF, sizeof (TurbConstVecs) );
    if (!turb_vecs) goto no_memory;


    //end allocation

    //Initialise turb constants
    char turb_pos[40];
    sprintf(turb_pos, "%s/turbulence.pos", work_dir);


    if (TurbConstsLoad 
        (   
			K,
            turb_pos,
            turb,
            turb_vecs
        )
	   ) { goto free_memory; }
//TODO add exceptions in function
    //TODO check headers

    //initialise turb constants ends

    //file reading starts
    fprintf(stderr, "\nReading data files...     ");
    /*      loop over all files        */

    Tensor2 meanS = {{0}};
    for (unsigned t = 0, n = 0;  t < K.iteration_num; t += K.t_gap) 
    {
        assert (n < K.SnapshotNum);
        fprintf(stderr, "\n%lf%%", 100.0*(double)n/(double)K.SnapshotNum); 
        fflush(NULL);
        //TODO load turbulence.pos only once!!

        char fname[60];
        sprintf (fname, "%s/t-%d.pos", work_dir, t);

       if ( MDLoadPos
            (
             K,
             fname,
             positions[n] 
            )
           ) { goto free_memory; }

       sprintf (fname, "%s/t-%d.pos", work_dir, t);
       if ( MDLoadDir
            (
             K,
             fname,
             positions[n] 
            )
           ) { goto free_memory; }
       
       sprintf (fname, "%s/turbField-%d.pos", work_dir, t);
       if ( MDLoadTurb
            (
             K,
             fname,
             turb_velocities[n] 
            )
           ) { goto free_memory; }


        if ( InitializeTurbModes
             (
                K,
                positions[n],
                turb_vecs,
                turb,
                kraich_modes[n],
                t
             )
           ) { goto free_memory; }
/*        if ( n > 0)
        {
            if ( InitializeTurbVelocities
                 (
                   K,
                   turb_vecs,
                   kraich_modes[n],
                   turb_velocities[n] 
                 )
               ) { goto free_memory; }
        }
*/
//FIXME generated turb velocities are not the same as the ones contained in the .pos!

//init strain rate tensor
        for (unsigned i = 0; i < K.PartNum; ++i)
        {
            InitializeStrainRateTensor
                (
                     strain_rate[n][i],
                     K,
                     turb_vecs,
                     kraich_modes[n][i]
                );
        } 
          n++;
    }
    //file reading ends
    //calc kinetic energy
    double mean_kinetic_energy =
        MeanKineticEnergy 
        (
         (const MDConstants) K,
         (const Molecule**) positions,
         (const TurbField**) turb_velocities
        );

    if( MeanStrainRateTensor
        (
             (const Tensor2**) strain_rate,
             (const MDConstants) K,
             meanS
        )
      ) goto free_memory;

    for (unsigned k = 0; k < kDIM; ++k)
    {
        char message[50] = {0};
        sprintf(message, "S[%u][%%u] = %%f\t", k);
        PrintVals (message, meanS[k], kDIM);
    }

    fprintf (stdout, "Mean kinetic energy = %lf", mean_kinetic_energy);
	goto free_memory;

//EXCEPTIONS
	no_memory:
	{
		fprintf(stderr, "no memory\n");
		goto free_memory;
	}

    free_memory:
	{
		for (unsigned n = 0; n < K.SnapshotNum; n++) 
		{
			free (positions[n]);
			positions[n] = NULL;

			free (turb_velocities[n]);
			turb_velocities[n] = NULL;

			free (strain_rate[n]);
			strain_rate[n] = NULL;


			for (unsigned i = 0; i < K.PartNum; ++i)
			{
				free (kraich_modes[n][i]);
				kraich_modes[n][i] = NULL;
			}

			free (kraich_modes[n]);
			kraich_modes[n] = NULL;
		}

		free(positions);
		positions = NULL;	

		free(turb_velocities);
		turb_velocities = NULL;	

		free(strain_rate);
		strain_rate = NULL;	

		free (turb);
		turb = NULL;

		free (turb_vecs);
		turb_vecs = NULL;

		free (kraich_modes);
		kraich_modes = NULL;

		goto exit;
	}

    exit:
    {
        return 0;
    }
	error:
	{
		return -1;
	}
}

