#ifndef MDLOAD_H
#define MDLOAD_H

#include "Molecule.h" //Molecule
#include "Turbulence.h" //TurbConsts, TurbConstVecs






extern char work_dir[];

int
TurbConstsLoad 
(
    const MDConstants K,
    TurbConsts* turb,
    TurbConstVecs* turb_vecs
);

int
MDLoadDir
(
    const MDConstants K, 
    Molecule* molecule, 
    const unsigned next_loop
);


int
MDLoadPos
(
    const MDConstants K, 
    Molecule* molecule, 
    const unsigned next_loop
);

int
MDLoadTurb
(
    const MDConstants K, 
    TurbField* turb_field,
    const unsigned next_loop
);

int
MDLoadMol
(
    const MDConstants K, 
	const char fname[],
	const unsigned pre_offset,
	const unsigned post_offset,
    Molecule* molecule
);

#endif /* MDLOAD_H */
