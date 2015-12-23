#include <assert.h> //assert()
#include <math.h> // pow()
#include <stdio.h> //fprintf()

#include "MDPostProcessing.h"

double
MeanKineticEnergy
(
    MDConstants K,
    Molecule **positions,
	TurbField **turb_velocities
    )
{
    double time_mean = 0;

    for (unsigned n = 0; n < K.SnapshotNum; ++n) // calc mean over time
    {
        double part_mean = 0; //for every new timestep

        for (unsigned i = 0; i < K.PartNum; ++i)//calc mean over particles
        {		
            double v[kDIM]; //temp array
            double v_squared = 0;

            for (unsigned j = 0; j < kDIM; ++j)
            {
                v[j] = K.v_0 * positions[n][i].direction[j] + turb_velocities[n][i].direction[j];
                v_squared += pow (v[j], 2);
            }

            double kin_energy = v_squared / 2.0; // mass is assumed = 1

            part_mean += kin_energy;	assert (part_mean >= 0);
        }

        part_mean /= K.PartNum;		

        time_mean += part_mean; assert (time_mean > 0);
    }	

    time_mean /= (K.iteration_num * K.delta_t);

    return time_mean;
}

