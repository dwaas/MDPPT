#ifndef TURBULENCE_H
#define TURBULENCE_H

//1 for each mode, absolute values
struct turb_var 
{
  double k;
  double omega;
  double delk;
  double cabs;
};

typedef struct turb_var TurbConsts;

//only for kraichnan's method.
//1 for each mode
struct turb_vecs_var
{
    double kn[kDIM]; // k_n unit vectors
    double c1n[kDIM]; //a_n cross k_n
    double c2n[kDIM]; //b_n cross k_n
    double d1[kDIM]; //c_1 cross k_n
    double d2[kDIM]; //c_2 cross k_n
    //TODO consistency in names
};
typedef struct turb_vecs_var TurbConstVecs;


struct kraichnan_mode_var
{
    double sin;
    double cos;
};
typedef struct kraichnan_mode_var KraichnanMode;


struct turb_field_var
{
    double direction[kDIM]; // velocity components of the field
    double vorticity[kDIM]; //vorticity components of the field
};

typedef struct turb_field_var TurbField;

typedef double Tensor2[kDIM][kDIM];

int 
InitializeTurbModes
(
     MDConstants K,
     Molecule* molecule,
     TurbConstVecs* turb_vecs,
     TurbConsts* turb,
     KraichnanMode** kraich_modes,
     unsigned itime
);

int
InitializeTurbVelocities
(
	MDConstants K,
	TurbConstVecs* turb_vecs,
	KraichnanMode** kraich_modes,
	TurbField* turb_velocities
);

//TODO call it initialize
void
StrainRateTensor
(
	Tensor2 S,
	const MDConstants K,
	const TurbConstVecs* turb_vecs,
	const KraichnanMode* modes
);

void 
MeanStrainRateTensor
(
	Tensor2** S,
	MDConstants K,
    Tensor2 meanS	
);


//helper functions
double 
DotProd 
(
	const double vec1[],
	const double vec2[]
);

void
NormalizeVector
(	
	const double consts[],
	const double in_vec[],
	double out_vec[]
);


#endif /* TURBULENCE_H */
