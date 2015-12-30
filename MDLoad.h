#ifndef MDLOAD_H
#define MDLOAD_H

#include "Molecule.h" //Molecule
#include "Turbulence.h" //TurbConsts, TurbConstVecs







int
TurbConstsLoad 
(
    const MDConstants K,
    const char fname[],
    TurbConsts* turb,
    TurbConstVecs* turb_vecs
);

int
MDLoadDir
(
    const MDConstants K, 
    const char fname[],
    Molecule* molecule 
);


int
MDLoadPos
(
    const MDConstants K, 
    const char fname[],
    Molecule* molecule 
);

int
MDLoadTurb
(
    const MDConstants K, 
    const char fname[],
    TurbField* turb_field
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
