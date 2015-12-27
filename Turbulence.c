#include <assert.h> //assert()
#include <math.h> //cos (); sin ()
#include <stdlib.h> //free, calloc





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
			double pos_vec[kDIM] = {0};

			NormalizeVector (K.L, molecule[i].position, pos_vec);
			double kn_dot_x = DotProd ( turb_vecs[modeIndex].kn, pos_vec);

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

 
void
StrainRateTensor
(
	Tensor2 S,
	const MDConstants K,
	const TurbConstVecs* turb_vecs,
	const KraichnanMode* modes
)
{
	for (unsigned f = 0; f < K.NF; ++f)
	{
		for (unsigned i = 0; i < kDIM; ++i)
		{
			for (unsigned j = 0; j < kDIM; ++j)
			{
				S[i][j] -= turb_vecs->c1n[i] * turb_vecs->kn[j] * modes[f].sin;
				S[i][j] += turb_vecs->c2n[i] * turb_vecs->kn[j] * modes[f].cos;
			}
		}
	}
	return;
}

void 
MeanStrainRateTensor
(
	Tensor2** S,
	MDConstants K,
    Tensor2 meanS	
)
{
	Tensor2* tempS = (Tensor2*) calloc (K.SnapshotNum, sizeof (Tensor2) );
	if (!tempS) goto no_memory;


	for (unsigned t = 0; t < K.SnapshotNum; ++t)
	{
		for (unsigned k = 0; k < K.PartNum; ++k)
		{
			// sum tensor
			for (unsigned i = 0; i < kDIM; ++i)
			{
				for (unsigned j = 0; j < kDIM; ++j)
				{
					tempS[t][i][j] += S[t][k][i][j];
				}
			}
		}
		//divide tensor
		for (unsigned i = 0; i < kDIM; ++i)
		{
			for (unsigned j = 0; j < kDIM; ++j)
			{
				tempS[t][i][j] /= K.PartNum;
			}
		}

		// sum tensor
		for (unsigned i = 0; i < kDIM; ++i)
		{
			for (unsigned j = 0; j < kDIM; ++j)
			{
				meanS[i][j] += tempS[t][i][j];
			}
		}

		//divide tensor
		for (unsigned i = 0; i < kDIM; ++i)
		{
			for (unsigned j = 0; j < kDIM; ++j)
			{
				meanS[i][j] /= K.SnapshotNum;
			}
		}
	}

free_memory:
	{
		free(tempS);
		tempS = NULL;	
	}

	//TODO assert that meanS is initialized to 0
	return;

	//EXCEPTIONS
no_memory:
	{
		fprintf(stderr, "no memory\n");
		goto free_memory;
	}
}

int
InitializeTurbVelocities
(
	MDConstants K,
	TurbConstVecs* turb_vecs,
	KraichnanMode** kraich_modes,
	TurbField* turb_velocities
)
{
	for (unsigned i = 0; i < K.PartNum; ++i)
	{
		double vel[kDIM] = {0};//temp array
		for (unsigned j = 0; j < kDIM; ++j)
		{
			for (unsigned f = 0; f < K.NF; ++f)
			{
				vel[j] += turb_vecs[f].c1n[j] * kraich_modes[i][f].cos;
				vel[j] += turb_vecs[f].c2n[j] * kraich_modes[i][f].sin;
			}	
			if (turb_velocities[i].direction[j] != vel[j])
			{
				goto error;
			}
		}	
	}
	return 0;
//TODO make this a test case
	error:
	{
		fprintf 
		(
			stderr,
			"the algorithm doesn't generate the correct velocities"
		);
		return -1;
	}
}


double 
DotProd 
(
	const double vec1[],
	const double vec2[]
)
{
	double dot_prod = 0;
	for (unsigned j = 0; j < kDIM; ++j)
	{
		dot_prod += vec1[j] * vec2[j];
	}//positions normalized on box length

	return dot_prod;
}

void
NormalizeVector
(	
	const double consts[],
	const double in_vec[],
	double out_vec[]
)
{
	for (unsigned j = 0; j < kDIM; ++j)
	{
		out_vec[j] = in_vec[j] * consts[j];
	}

	return;
}


