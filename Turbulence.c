#include <assert.h> //assert()
#include <math.h> //cos (); sin ()





#include "MDConstants.h" //MDConstants
#include "Molecule.h" //Molecule
#include "Turbulence.h" //KraichnanMode

//TODO MDInitialize vorticities
int
InitializeTurbModes
(
	MDConstants K,
	Molecule* molecule,
	TurbConstVecs* turb_vecs,
	TurbConsts* turb,
	KraichnanMode** kraich_modes,
	unsigned t
)
{
	assert (t < K.iteration_num);

	//TODO profile & parallelise
	//TODO check valid input/output 
/*
#pragma omp parallel for default (none) \
shared (K.PartNum, molecule, turb_vecs, turb, K.Lx, K.Ly, K.Lz, t, delta_t, NF)\
schedule(dynamic, 5000)
*/
	for(unsigned i = 0; i < K.PartNum; ++i)
	{
		for(unsigned modeIndex = 0; modeIndex < K.NF; ++modeIndex)
		{
			double kn_dot_x = 0;

			for (unsigned j = 0; j < kDIM; ++j)
			{
				kn_dot_x += turb_vecs[modeIndex].kn[j] * molecule[i].position[j] * K.L[j];
			}//positions normalized on box length

			unsigned real_time = t * K.delta_t;                

			//in LaTeX would be:  \Omega_n = |\kappa| (\vec k \cdot \vec x) + |\omega| t
			double Omega_n = turb[modeIndex].k * kn_dot_x +
				turb[modeIndex].omega * real_time;

			//cos and sin
			kraich_modes[i][modeIndex].sin = sin (Omega_n);
			kraich_modes[i][modeIndex].cos = cos (Omega_n);
		}
	}
	return 0;
}
