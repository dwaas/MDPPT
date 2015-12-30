#include <assert.h> //assert
#include <stdio.h> // FILE; fprintf(), fread()


#include "debug.h"

#include "MDConstants.h" //MDConstants
#include "MDLoad.h" 
#include "Molecule.h" // Molecule
#include "Turbulence.h" // TurbConsts; TurbConstVecs;



//TODO const pointers vs copies
int
TurbConstsLoad 
(
    const MDConstants K,
    const char fname[],
 	TurbConsts* turb,
 	TurbConstVecs* turb_vecs
)
{
   	FILE* fp;
    check
	(
		 (fp = fopen (fname, "r")),
		 "\nError opening file %s\n", 
		 fname
	);

    unsigned count_scan = 0;

    for(unsigned i = 0; i < K.NF; ++i)
	{
		count_scan += fread(&turb[i].k, sizeof(double), 1, fp);
		count_scan += fread(&turb[i].omega,  sizeof(double), 1, fp);
		count_scan += fread(&turb[i].delk, sizeof(double), 1, fp);
		count_scan += fread(&turb[i].cabs, sizeof(double), 1, fp);
		count_scan += fread (&turb_vecs[i].kn, sizeof(double), kDIM, fp);
		count_scan += fread (&turb_vecs[i].c1n, sizeof(double), kDIM, fp);
		count_scan += fread (&turb_vecs[i].c2n, sizeof(double), kDIM, fp);
//TODO fread vs fscanf
	}
	
	check 
	(
		(count_scan == K.NF * 13),
        "Wrong number of turbulence points in:\
		\n%s\
		\nnum = %u",
		fname,
		count_scan
	);
    //TODO macro to get the third argument of fread

   	fprintf (stderr, "\n%s read.\n", fname);
    fclose(fp);

	return 0;

	error:
	{
		fclose(fp);
		return -1;
    }
}



int
MDLoadPos
(
    const MDConstants K, 
    const char fname[],
    Molecule* molecule
)
{
	//TODO check that K is initialised; check input is valid.

   	FILE* fp;
    check
	(
		 (fp = fopen (fname, "r")),
		 "\nError opening file %s\n", 
		 fname
	);
		
    unsigned count_scan = 0;
//TODO find tests for directions and turb velocities too

    double pos[kDIM];//temp variable
	double dump[kDIM];
    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        //TODO disk I/O checks

		count_scan += fread(&pos, sizeof(double), kDIM, fp);
        for (unsigned j = 0; j < kDIM; ++j)
        {
			check 
			(
				!(pos[j] < K.L[j] / (-2.0) || pos[j] > K.L[j] / 2.0),
				"The particles in the %u-th position\
                are not in the expected range for this simulation,\
                please check\
				\n%s\n", 
				j, 
				fname
			);

            molecule[i].position[j] = K.side_minus1[j] * pos[j];
        }// molecule.position[] and molecule.direction[] are initialized		
		fread(&dump, sizeof(double), kDIM, fp);
    } 



	check 
	(
		(count_scan == K.PartNum * kDIM),
        "Wrong number of turbulence points in:\
		\n%s\
		\nnum = %u",
		fname,
		count_scan
	);
 
	fprintf (stderr, "\n%s read.\n", fname);
    fclose(fp);
	
	return 0;

 	error:
	{
		fclose(fp);
		return -1;
    }      
}

int
MDLoadTurb 
(
    const MDConstants K, 
    const char fname[],
	TurbField* turb_field
)
{

   	FILE* fp;
    check
	(
		 (fp = fopen (fname, "r")),
		 "\nError opening file %s\n", 
		 fname
	);
		
	
    unsigned count_scan = 0;
    for (unsigned i = 0; i < K.PartNum; ++i)
    {
		count_scan += fread(&turb_field[i], sizeof(double), kDIM, fp);
    } 


	check 
	(
		(count_scan == K.PartNum * kDIM),
        "Wrong number of turbulence points in:\
		\n%s\
		\nnum = %u",
		fname,
		count_scan
	);
 


	fprintf (stderr, "\n%s read.\n", fname);
    fclose(fp);

	return 0;

 	error:
	{
		fclose(fp);
		return -1;
    }    
}

int
MDLoadDir
(
    const MDConstants K, 
    const char fname[],
    Molecule* molecule
)
{
	//TODO check that K is initialised; check input is valid.
   	FILE* fp;
    check
	(
		 (fp = fopen (fname, "r")),
		 "\nError opening file %s\n", 
		 fname
	);
		
    unsigned count_scan = 0;
//TODO find tests for directions and turb velocities too

	double dump[kDIM];
    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        //TODO disk I/O checks

		fread(&dump, sizeof(double), kDIM, fp);
		count_scan += fread(&molecule[i].direction, sizeof(double), kDIM, fp);
    } 

	check 
	(
		(count_scan == K.PartNum * kDIM),
        "Wrong number of turbulence points in:\
		\n%s\
		\nnum = %u",
		fname,
		count_scan
	);
 
	fprintf (stderr, "\n%s read.\n", fname);
    fclose(fp);
	
	return 0;

	error:
	{
		fclose(fp);
		return -1;
    }
}
//FIXME
//TODO refactor and extend
int
MDLoadMol
(
    const MDConstants K, 
	const char fname[],
	const unsigned pre_offset,
	const unsigned post_offset,
    Molecule* molecule
)
{
	assert (pre_offset < 2);
	assert (post_offset < 2);
	//TODO check that K is initialised; check input is valid.
   	FILE* fp;
    check
	(
		 !(fp = fopen (fname, "r")),
		 "\nError opening file %s\n", 
		 fname
	);
	
    unsigned count_scan = 0;
//TODO find tests for directions and turb velocities too

	double dump[kDIM];
    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        //TODO disk I/O checks

        fread(&dump, sizeof(double), pre_offset*kDIM, fp);
		count_scan += fread(&molecule[i].position, sizeof(double), kDIM, fp);
        count_scan += fread(&dump, sizeof(double), post_offset*kDIM, fp);
		} 


    fclose(fp);

	check 
	(
		!(count_scan == K.PartNum * kDIM),
        "Wrong number of turbulence points in:\
		\n%s\
		\nnum = %u",
		fname,
		count_scan
	);
 
	fprintf (stderr, "\n%s read.\n", fname);
	
	return 0;

	error:
	{
		fclose(fp);
		return -1;
    }

}



