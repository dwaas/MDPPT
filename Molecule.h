#ifndef MOLECULE_H
#define MOLECULE_H

struct molecule_var
{
	double x,y,z;
	double torque_x,torque_y,torque_z;
	double e_x,e_y,e_z;
	double lambda;
	unsigned neighbors;
};

typedef struct molecule_var Molecule;

#endif /* MOLECULE_H */
