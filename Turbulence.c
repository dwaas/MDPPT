#include <assert.h> //assert()
#include <math.h> //cos (); sin ()
#include <stdlib.h> //free, calloc





#include "Allocate.h" //FREE
#include "debug.h"
#include "MDConstants.h" //MDConstants
#include "Molecule.h" //Molecule
#include "Turbulence.h" //KraichnanMode

//TODO MDInitialize vorticities
int
InitializeTurbModes
(
	const MDConstants K,
	const Molecule* molecule,
	const TurbConstVecs* turb_vecs,
	const TurbConsts* turb,
	KraichnanMode** kraich_modes,
	const unsigned t
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
InitializeStrainRateTensor
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

int 
MeanStrainRateTensor
(
	const Tensor2** S,
	const MDConstants K,
    Tensor2 meanS	
)
{
	Tensor2* tempS = NULL;
    ALLOC (tempS, K.SnapshotNum); 

	for (unsigned t = 0; t < K.SnapshotNum; ++t)
	{
		MeanTensor ( tempS[t], S[t], K.PartNum );
	}
	MeanTensor ( meanS, tempS, K.SnapshotNum );


	//TODO assert that meanS is initialized to 0

    FREE (tempS);
	return 0;

	error:
	{
		FREE (tempS);
		return -1;
	}
}

int
InitializeTurbVelocities
(
	const MDConstants K,
	const TurbConstVecs* turb_vecs,
	const KraichnanMode** kraich_modes,
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


void
SumTensor
(
	Tensor2 tempS,
	const Tensor2 S
)
{
	for (unsigned i = 0; i < kDIM; ++i)
	{
		for (unsigned j = 0; j < kDIM; ++j)
		{
			tempS[i][j] += S[i][j];
		}
	}
	return;
}
void
DivideTensor
(
	Tensor2 S,
	const double K
)
{
	//divide tensor
	for (unsigned i = 0; i < kDIM; ++i)
	{
		for (unsigned j = 0; j < kDIM; ++j)
		{
			S[i][j] /= K;
		}
	}
	return;
}

void
MeanTensor
(	
	Tensor2 meanS,
	const Tensor2* S,
	const unsigned size
)
{
	for (unsigned k = 0; k < size; ++k)
	{
		SumTensor (meanS, S[k]);
	}
	DivideTensor (meanS, size);

	return;
}

