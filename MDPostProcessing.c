#include <assert.h> //assert()
#include <math.h> // pow()
#include <stdlib.h> //exit()
#include <stdio.h> //fprintf()

#include "MDPostProcessing.h"

double
MeanKineticEnergy
                (
                    MDConstants K,
                    Molecule **positions
                )
{
	double time_mean = 0;
	for (unsigned n = 1; n < K.SnapshotNum; ++n) // calc mean over time
	{
		double part_mean	= 0; //for every new timestep
		
		for (unsigned i = 0; i < K.PartNum; ++i)//calc mean over particles
		{
                double v_squared = 
                                    pow (positions[n][i].e_x, 2) +
                                    pow (positions[n][i].e_y, 2) +
                                    pow (positions[n][i].e_z, 2); 
 
                double kin_energy = v_squared / 2.0; // mass is assumed = 1

				part_mean += kin_energy;	assert (part_mean >= 0);
		}

		part_mean /= K.PartNum;		

		time_mean += part_mean; assert (time_mean > 0);
	}	
	
    time_mean /= (K.iteration_num * K.delta_t);
    time_mean *= K.v_0; // modulus of each velocity

    return time_mean;
}
double
MeanStrainRateTensor
                    (
                        MDConstants K,
                        Molecule** positions
                    )
{
	double time_mean = 0;
	for (unsigned n = 1; n < K.SnapshotNum; ++n) // calc mean over time
	{
		double part_mean	= 0; //for every new timestep
		
		for (unsigned i = 0; i < K.PartNum; ++i)//calc mean over particles
		{
                double v_squared = 
                                    pow (positions[n][i].e_x, 2) +
                                    pow (positions[n][i].e_y, 2) +
                                    pow (positions[n][i].e_z, 2); 
 
                double kin_energy = v_squared / 2.0; // mass is assumed = 1

				part_mean += kin_energy;	assert (part_mean >= 0);
		}

		part_mean /= K.PartNum;		

		time_mean += part_mean; assert (time_mean > 0);
	}	
	
    time_mean /= (K.iteration_num * K.delta_t);
    time_mean *= K.v_0; // modulus of each velocity
    return 1.0;
//     return MeanStrainRateTensor;
}

