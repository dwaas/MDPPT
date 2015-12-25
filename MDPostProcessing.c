#include <assert.h> //assert()
#include <math.h> // pow()
#include <stdio.h> //fprintf()

#include "MDPostProcessing.h"

double
FunctionAverage
(
    const double array[],
    const unsigned size,
    double (*f) (double)
);

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
            double v_part[kDIM] = {0};
            double v_0[kDIM];  
            InitConstArray (v_0, kDIM, K.v_0);
			// v_0 * \vec{e_part}
            NormalizeVector (v_0, positions[n][i].direction, v_part);
            double v[kDIM]; //temp array
            SumVector (v_part, turb_velocities[n][i].direction, v);
             
			double kin_energy = KineticEnergy (v);

            part_mean += kin_energy;	assert (part_mean >= 0);
        }

        part_mean /= K.PartNum;		

        time_mean += part_mean; assert (time_mean > 0);
    }	

    time_mean /= (K.iteration_num * K.delta_t);

    return time_mean;
}

//TODO use fma()

double
FunctionAverage
(
    const double array[],
    const unsigned size,
    double (*f) (double)
)
{   
    double sum = 0;   
    for (unsigned i = 0; i < size; ++i)
    {
        sum += f(array[i]);

    }

    return sum / size;
}

double
KineticEnergy
(
    const double v[kDIM]
)
{
    double v_squared = 0;
    for (unsigned j = 0; j < kDIM; ++j)
    {
        v_squared += pow (v[j], 2);
    }

    return v_squared / 2.0; // MASS is assumed to be 1
}

//TODO limit size of in_vecs
void
SumVector
(
    const double in_vec1[],
    const double in_vec2[],
    double out_vec[]
)
{
    for (unsigned j = 0; j < kDIM; ++j)
	{
		out_vec[j] = in_vec1[j] + in_vec2[j];
	}
    return;
}

void
InitConstArray
(
    double vec[],
    unsigned size,
    const double k
)
{
    for (unsigned j = 0; j < size; ++j)
    {
        vec[j] = k;
    }
    return;
}
