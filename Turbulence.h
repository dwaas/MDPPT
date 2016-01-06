#ifndef TURBULENCE_H
#define TURBULENCE_H
//TODO avoid to forget absolute values
//1 for each mode, absolute values
struct turb_var 
{
  double k;
  double omega;
  double delk;
  double cabs;
};

typedef struct turb_var TurbConsts;

//1 for each mode
//all turb_vecs are unit vectors (and cross products of them
//they need to be multiplied by their modulus found in turb_var
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
     const MDConstants K,
     const Molecule* molecule,
     const TurbConstVecs* turb_vecs,
     const TurbConsts* turb,
     KraichnanMode** kraich_modes,
     const unsigned itime
);

int
InitializeTurbVelocities
(
	const MDConstants K,
	const TurbConstVecs* turb_vecs,
	const TurbConsts* turb,
	const KraichnanMode** kraich_modes,
	TurbField* turb_velocities
);

//TODO call it initialize
void
InitializeStrainRateTensor
(
	Tensor2 S,
	const MDConstants K,
	const TurbConstVecs* turb_vecs,
	const TurbConsts* turb,
	const KraichnanMode* modes
);

int 
MeanStrainRateTensor
(
	const Tensor2** S,
	const MDConstants K,
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

void
DivideTensor
(
	Tensor2 S,
	const double K
);
void
SumTensor
(
	Tensor2 tempS,
	const Tensor2 S
);


void
MeanTensor
(	
	Tensor2 meanS,
	const Tensor2* S,
	const unsigned size
);


#endif /* TURBULENCE_H */
