#include <assert.h> //assert()
#include <math.h> //cos (); sin ()





#include "MDConstants.h" //MDConstants
#include "Molecule.h" //Molecule
#include "Turbulence.h" //KraichnanMode

//TODO MDInitialize vorticities
void
InitializeTurbModes
					(
						MDConstants K,
						Molecule* molecule,
    					TurbConstVecs* turb_vecs,
    					TurbConsts* turb,
						KraichnanMode** kraich_modes,
						unsigned itime)
{
	assert (itime < K.iteration_num);

	//TODO profile & parallelise
	//TODO check valid input/output 
/*
	#pragma omp parallel for default (none) \
	shared (K.PartNum, molecule, turb_vecs, turb, K.Lx, K.Ly, K.Lz, itime, delta_t, NF)\
		schedule(dynamic, 5000)
*/
		for(unsigned i = 0; i < K.PartNum; ++i)
		{
//TODO/R make more readable
			double 
				xPart = molecule[i].x,
				yPart = molecule[i].y,
				zPart = molecule[i].z;

			for(unsigned modeIndex = 0; modeIndex < K.NF; ++modeIndex)
			{
				double kn_x = turb_vecs[modeIndex].kn_x,
					   kn_y = turb_vecs[modeIndex].kn_y,
					   kn_z = turb_vecs[modeIndex].kn_z;

				// \Omega_n = |\kappa| (\vec k \cdot \vec x) + |\omega| t
				double Omega_n = turb[modeIndex].k * 
								(kn_x * xPart * K.Lx + kn_y * yPart * K.Lx +  kn_z * zPart * K.Lx) +
								turb[modeIndex].omega * 
								(itime-1) * K.delta_t; 
				//cos and sin
				kraich_modes[i][modeIndex].sin = sin (Omega_n);
				kraich_modes[i][modeIndex].cos = cos (Omega_n);

			}
		}
}
