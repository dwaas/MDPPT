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
                    MDConstants K,
                    Molecule** positions,
                    TurbField** turb_velocities
                );
double
MeanStrainRateTensor
                    (
                        MDConstants K,
                        Molecule** positions
                    );
#endif /* MDPOSTPROCESSING_H */
