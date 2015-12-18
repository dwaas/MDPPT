#ifndef TURBULENCE_H
#define TURBULENCE_H

//1 for each mode
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
  double kn_x, kn_y, kn_z; // k_n unit vectors
  double c1n_x, c1n_y, c1n_z; //a_n cross k_n
  double c2n_x, c2n_y, c2n_z; //b_n cross k_n
  double d1_x, d1_y, d1_z; //c_1 cross k_n
  double d2_x, d2_y, d2_z; //c_2 cross k_n
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
  double e_x, e_y, e_z; // velocity components of the field
  double vort_x, vort_y, vort_z; //vorticity components of the field
};

typedef struct turb_field_var TurbField;

void 
InitializeTurbModes
					(
						MDConstants K,
						Molecule* molecule,
    					TurbConstVecs* turb_vecs,
    					TurbConsts* turb,
						KraichnanMode** kraich_modes,
						unsigned itime
                    );
#endif /* TURBULENCE_H */
