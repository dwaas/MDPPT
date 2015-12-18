#include <assert.h> //assert()
#include <stdbool.h> //true, false
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>



#include "MDConstants.h" //MDConstants;
#include "MDLoad.h" // MDrestart();
#include "MDPostProcessing.h" // MeanKineticEnergy()
#include "Molecule.h" // Molecule
#include "Turbulence.h" //TurbConsts; TurbConstVecs

#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))

//TODO documentation



char work_dir[40]; //TODO find a safer way

int
main(int argc, char *argv[])
{
    Molecule** positions;
    TurbField** turb_velocities;
    TurbConsts* turb;
    TurbConstVecs* turb_vecs;


    //ARGUMENT CHECK	
    if (argc != 2) //TODO improve interface 
    {
        fprintf(stderr, "\nError! Number of arguments is wrong! usage: ./post_proc working_directory_path!!!!\n\n");
        exit(1);
    }
    //CONSTANTS
    MDConstants K;
    MDConstants* K_ptr = &K;
    Initialize (K_ptr, argv[1]);//work_dir is now initialised

    //allocate 
    positions = (Molecule**) calloc (K.SnapshotNum, sizeof (Molecule*) );
    if (!positions) fprintf(stderr, "no memory\n");

    turb_velocities = (TurbField**) calloc (K.SnapshotNum, sizeof (TurbField*) );
    if (!turb_velocities) fprintf(stderr, "no memory\n");

    KraichnanMode*** kraich_modes = (KraichnanMode***) calloc (K.SnapshotNum, sizeof (KraichnanMode**) );
    if (!kraich_modes) fprintf (stderr, "no memory\n");

    for (unsigned n = 0; n < K.SnapshotNum; ++n) 
    {
        positions[n] = (Molecule*) calloc (K.PartNum, sizeof (Molecule) );
        if (!positions[n]) fprintf (stderr, "no memory\n");

        turb_velocities[n] = (TurbField*) calloc (K.PartNum, sizeof (TurbField) );
        if (!turb_velocities[n]) fprintf (stderr, "no memory\n");

        kraich_modes[n] = (KraichnanMode**) calloc (K.PartNum, sizeof (KraichnanMode*) );
        if (!kraich_modes[n]) fprintf (stderr, "no memory\n");


        for (unsigned i = 0; i < K.PartNum; ++i)
        {
            kraich_modes[n][i] = (KraichnanMode*) calloc (K.NF, sizeof (KraichnanMode) );
            if (!kraich_modes[n][i]) fprintf (stderr, "no memory\n");

        }

    }

    turb = (TurbConsts*) calloc (K.NF, sizeof (TurbConsts) );
    if (!turb) fprintf(stderr, "no memory\n");

    turb_vecs = (TurbConstVecs*) calloc (K.NF, sizeof (TurbConstVecs) );
    if (!turb_vecs) fprintf(stderr, "no memory\n");


    //end allocation

    //Initialise turb constants
    TurbConstsLoad (
            K,
            turb,
            turb_vecs
            );
    //TODO check headers

    //initialise turb constants ends

    //file reading starts
    fprintf(stderr, "\nReading data files...     ");
    /*      loop over all files        */

    for (unsigned t = 0, n = 0;  t < K.iteration_num; t += K.t_gap) 
    {
        assert (n < K.SnapshotNum);
        fprintf(stderr, "\n%lf%%", 100.0*(double)n/(double)K.SnapshotNum); 
        fflush(NULL);
        //TODO load turbulence.pos only once!!

        MDLoad 
            (
             K,
             positions[n],
             turb_velocities[n], 
             t 
            );

        InitializeTurbModes
            (
                K,
                positions[n],
                turb_vecs,
                turb,
                kraich_modes[n],
                t
             );
        n++;
    }
    //file reading ends
    //calc kinetic energy
    double mean_kinetic_energy =
        MeanKineticEnergy 
        (
         K,
         positions
        );

    fprintf (stdout, "Mean kinetic energy = %lf", mean_kinetic_energy);
    //free memory 
    for (unsigned n = 0; n < K.SnapshotNum; n++) 
    {
        free (positions[n]);
        positions[n] = NULL;

        free (turb_velocities[n]);
        turb_velocities[n] = NULL;

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

    free (turb);
    turb = NULL;

    free (turb_vecs);
    turb_vecs = NULL;

    free (kraich_modes);
    kraich_modes = NULL;
    //free memory ends
    return 0;

}
