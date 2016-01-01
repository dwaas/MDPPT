#include <assert.h> //assert()
#include <stdio.h>

#include "Allocate.h" //FREE
#include "debug.h" //check(), check_mem()
#include "MDConstants.h" //MDConstants;
#include "MDLoad.h" // MDLoad();
#include "MDPostProcessing.h" // MeanKineticEnergy()
#include "Molecule.h" // Molecule
#include "Turbulence.h" //TurbConsts; TurbConstVecs


//TODO documentation
//TODO testsuite: wrong args, wrong input file, wrong data in input file, no memory

int
main
(
	int argc, 
	char *argv[]
)
{
        //declare
        KraichnanMode*** kraich_modes = NULL;

        Molecule** positions = NULL;
        TurbField** turb_velocities = NULL;
        Tensor2** strain_rate = NULL;    

        TurbConsts* turb = NULL;
        TurbConstVecs* turb_vecs = NULL;


        check (
                        argc == 2,
                        "\nNumber of arguments is wrong! \
                        \nusage: ./post_proc working_directory_path!!!!\n\n"
              );

        char work_dir[40];
        sprintf (work_dir, "%s", argv[1]);

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
        ALLOC3D (kraich_modes, K.SnapshotNum, K.PartNum, K.NF);

        ALLOC2D (positions, K.SnapshotNum, K.PartNum);
        ALLOC2D (turb_velocities, K.SnapshotNum, K.PartNum);
        ALLOC2D (strain_rate, K.SnapshotNum, K.PartNum);

        ALLOC (turb, K.NF);
        ALLOC (turb_vecs, K.NF);

        //Initialise turb constants
        char turb_pos[40];
        sprintf(turb_pos, "%s/turbulence.pos", work_dir);

        check
                (
                 !TurbConstsLoad 
                 (   
                  K,
                  turb_pos,
                  turb,
                  turb_vecs
                 ),
                 "TurbConstLoad failed"
                );
        //TODO check headers

        //initialise turb constants ends

        //file reading starts
        fprintf(stderr, "\nReading data files...     ");

        Tensor2 meanS = {{0}};
        for (unsigned t = 0, n = 0;  t < K.iteration_num; t += K.t_gap) 
        {
                assert (n < K.SnapshotNum);
                fprintf(stderr, "\n%lf%%", 100.0*(double)n/(double)K.SnapshotNum); 
                fflush(NULL);
                //TODO load turbulence.pos only once!!

                char fname[60];
                sprintf (fname, "%s/t-%d.pos", work_dir, t);

                check
                        (
                         !MDLoadPos
                         (
                          K,
                          fname,
                          positions[n] 
                         ),
                         "MDLoadPos failed"
                        );

                sprintf (fname, "%s/t-%d.pos", work_dir, t);

                check
                        (
                         !MDLoadDir
                         (
                          K,
                          fname,
                          positions[n] 
                         ),
                         "MDLoadDir failed"
                        );
                sprintf (fname, "%s/turbField-%d.pos", work_dir, t);

                check
                        (
                         !MDLoadTurb
                         (
                          K,
                          fname,
                          turb_velocities[n] 
                         ),
                         "MDLoadTurb failed"
                        );

                check
                        (
                         !InitializeTurbModes
                         (
                          K,
                          positions[n],
                          turb_vecs,
                          turb,
                          kraich_modes[n],
                          t
                         ),
                         "InitializeTurbModes failed"
                        );
                /*
                   if ( n > 0)
                   {
                   check
                   ( 
                   !InitializeTurbVelocities
                   (
                   K,
                   turb_vecs,
                   (const KraichnanMode**) kraich_modes[n],
                   turb_velocities[n] 
                   ),
                   "InitializeTurbVelocities failed"
                   );
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

        check
                (

                 !MeanStrainRateTensor
                 (
                  (const Tensor2**) strain_rate,
                  (const MDConstants) K,
                  meanS
                 ),
                 "MeanStrainRateTensor failed"
                );

        for (unsigned k = 0; k < kDIM; ++k)
        {
                char message[50] = {0};
                sprintf(message, "S[%u][%%u] = %%f\t", k);
                PrintVals (message, meanS[k], kDIM);
        }

        fprintf (stdout, "Mean kinetic energy = %lf", mean_kinetic_energy);

        FREE3D (kraich_modes, K.PartNum, K.SnapshotNum);

        FREE2D (positions, K.SnapshotNum);
        FREE2D (turb_velocities, K.SnapshotNum);
        FREE2D (strain_rate, K.SnapshotNum);

        FREE (turb);
        FREE (turb_vecs);

        return 0;

error:
        {
                FREE3D (kraich_modes, K.PartNum, K.SnapshotNum);

                FREE2D (positions, K.SnapshotNum);
                FREE2D (turb_velocities, K.SnapshotNum);
                FREE2D (strain_rate, K.SnapshotNum);

                FREE (turb);
                FREE (turb_vecs);

                return -1;
        }
}

