#ifndef MDRESTART_H
#define MDRESTART_H

#include "Molecule.h" //Molecule
#include "Turbulence.h" //TurbConsts, TurbConstVecs






extern char work_dir[];

void MDLoad 
(
    const MDConstants K, 
    Molecule* molecule, 
    TurbField* turb_field,
    TurbConsts* turb,
    TurbConstVecs* turb_vecs,
    const unsigned next_loop
);

#endif /* MDRESTART_H */
