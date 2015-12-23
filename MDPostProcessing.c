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
		 //TODO refactor molecule with position ARRAY
				double v_x = K.v_0 * positions[n][i].e_x + turb_velocities[n][i].e_x;
				double v_y = K.v_0 * positions[n][i].e_y + turb_velocities[n][i].e_y;
				double v_z = K.v_0 * positions[n][i].e_z + turb_velocities[n][i].e_z;

                double v_squared = 
                                    pow (v_x, 2) +
                                    pow (v_y, 2) +
                                    pow (v_z, 2); 
 
                double kin_energy = v_squared / 2.0; // mass is assumed = 1

				part_mean += kin_energy;	assert (part_mean >= 0);
		}

		part_mean /= K.PartNum;		

		time_mean += part_mean; assert (time_mean > 0);
	}	
	
    time_mean /= (K.iteration_num * K.delta_t);

    return time_mean;
}
/*
double
MeanStrainRateTensor
                    (
                        MDConstants K,
                        Molecule** positions
                    )
{
	double time_mean = 0;
	for (unsigned n = 0; n < K.SnapshotNum; ++n) // calc mean over time
	{
		double part_mean	= 0; //for every new timestep
		
		for (unsigned i = 0; i < K.PartNum; ++i)//calc mean over particles
		{
			for (unsigned f = 0; f < K.NF; ++f)
			{
			
			}
		}

		part_mean /= K.PartNum;		

		time_mean += part_mean; assert (time_mean > 0);
	}	
	
    time_mean /= (K.iteration_num * K.delta_t);
    time_mean *= K.v_0; // modulus of each velocity
    return 0;
//     return MeanStrainRateTensor;
}
*/

