#ifndef MDPOSTPROCESSING_H
#define MDPOSTPROCESSING_H

#include "MDConstants.h"    //MDConstants
#include "Molecule.h"       //Molecule
#include "Turbulence.h"     //TurbField

//periodic boundary conditions
#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))

double
MeanKineticEnergy
(
    const MDConstants K,
    const Molecule** positions,
    const TurbField** turb_velocities
);

double
KineticEnergy
(
    const double v[kDIM]
);

void
SumVector
(
    const double in_vec1[],
    const double in_vec2[],
    double out_vec[]
);

void
InitConstArray
(
    double vec[],
    unsigned size,
    const double k
);

//TODO gsl_vector
#endif /* MDPOSTPROCESSING_H */
