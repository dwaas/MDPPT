#ifndef MOLECULE_H
#define MOLECULE_H

#include "MDConstants.h" //kDIM

struct molecule_var
{
    double position[kDIM];
    double torque[kDIM];
    double direction[kDIM];
    double lambda;
    unsigned neighbors;
};

typedef struct molecule_var Molecule;

#endif /* MOLECULE_H */
