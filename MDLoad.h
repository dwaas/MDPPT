#ifndef MDLOAD_H
#define MDLOAD_H

#include "Molecule.h" //Molecule
#include "Turbulence.h" //TurbConsts, TurbConstVecs






extern char work_dir[];

void
TurbConstsLoad (
			const MDConstants K,
    		TurbConsts* turb,
    		TurbConstVecs* turb_vecs
			);

void MDLoad 
(
    const MDConstants K, 
    Molecule* molecule, 
    TurbField* turb_field,
    const unsigned next_loop
);

#endif /* MDLOAD_H */
