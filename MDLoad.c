#include <assert.h> //assert
#include <stdio.h> // FILE; fprintf(), fread()

#include "MDConstants.h" //MDConstants
#include "MDLoad.h" 
#include "Molecule.h" // Molecule
#include "Turbulence.h" // TurbConsts; TurbConstVecs;


//u cross v = u12
void cross(double x1,double y1,double z1,double x2,double y2,double z2,double *x12,double *y12,double *z12){

  *x12 = y1*z2-y2*z1;
  *y12 = z1*x2-z2*x1;
  *z12 = x1*y2-x2*y1;

}


//TODO const pointers vs copies
int
TurbConstsLoad 
(
    const MDConstants K,
 	TurbConsts* turb,
 	TurbConstVecs* turb_vecs
)
{
	char fname[60];	
    sprintf(fname, "%s/turbulence.pos", work_dir);
   	FILE* fp = fopen(fname, "r");

    unsigned count_scan = 0;

    for(unsigned i = 0; i < K.NF; ++i)
	{
		count_scan += fread(&turb[i].k, sizeof(double), 1, fp);
		count_scan += fread(&turb[i].omega,  sizeof(double), 1, fp);
		count_scan += fread(&turb[i].delk, sizeof(double), 1, fp);
		count_scan += fread(&turb[i].cabs, sizeof(double), 1, fp);

		count_scan += fread(&turb_vecs[i].kn_x, sizeof(double), 1, fp);
		count_scan += fread(&turb_vecs[i].kn_y, sizeof(double), 1, fp);
		count_scan += fread(&turb_vecs[i].kn_z, sizeof(double), 1, fp);

		count_scan += fread(&turb_vecs[i].c1n_x, sizeof(double), 1, fp);
		count_scan += fread(&turb_vecs[i].c1n_y, sizeof(double), 1, fp);
		count_scan += fread(&turb_vecs[i].c1n_z, sizeof(double), 1, fp);

		count_scan += fread(&turb_vecs[i].c2n_x, sizeof(double), 1, fp);
		count_scan += fread(&turb_vecs[i].c2n_y, sizeof(double), 1, fp);
		count_scan += fread(&turb_vecs[i].c2n_z, sizeof(double), 1, fp);
/*
//calc d2
		cross (
					turb_vecs[i].c1n_x,
					turb_vecs[i].c1n_y,
					turb_vecs[i].c1n_z,
					turb_vecs[i].kn_x,
					turb_vecs[i].kn_y,
					turb_vecs[i].kn_z,
					&turb_vecs[i].d1_x,
					&turb_vecs[i].d1_y,
					&turb_vecs[i].d1_z
			);

	cross (
					turb_vecs[i].c2n_x,
					turb_vecs[i].c2n_y,
					turb_vecs[i].c2n_z,
					turb_vecs[i].kn_x,
					turb_vecs[i].kn_y,
					turb_vecs[i].kn_z,
					&turb_vecs[i].d2_x,
					&turb_vecs[i].d2_y,
					&turb_vecs[i].d2_z
			);
*/
//TODO refactor cross products!

	}
    fclose(fp);

    if(count_scan != K.NF * 13) 
    {
        fprintf (stderr, "Wrong number of turbulence points in %s.\n", fname);
		return -1;
    }
	fprintf (stderr, "\n%s read.\n", fname);

	return 0;
}



int
MDLoad 
(
    const MDConstants K, 
    Molecule* molecule, 
	TurbField* turb_field,
    const unsigned next_loop //current timestep
)
{
	//TODO check that K is initialised; check input is valid.
	assert (next_loop <= K.iteration_num);
//TODO check if positions are -L/2< x <+L/2

    double x1, y1, z1;

    char fname[60];
    FILE *fp;
    unsigned count_scan;
    double 
            side_x_minus1 = 1.0 / K.Lx,
            side_y_minus1 = 1.0 / K.Ly,
            side_z_minus1 = 1.0 / K.Lz;

	sprintf(fname,"%s/t-%d.pos",work_dir, next_loop);
    fp = fopen(fname, "r");
 	if ( !(fp = fopen(fname, "r")) )
		{
			fprintf(stderr, "\nError opening file %s;  Program aborted!\n\n", fname);
			return -1;
		}

   count_scan = 0;
//TODO find tests for directions and turb velocities too
    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        count_scan += fread(&x1, sizeof(double), 1, fp);
        count_scan += fread(&y1, sizeof(double), 1, fp);
        count_scan += fread(&z1, sizeof(double), 1, fp);
//TODO disk I/O checks
		
		if (x1 < K.Lx / (- 2.0) || x1 > K.Lx / 2.0) 
		{
			printf("\nThe particles x positions are not in the expexpected range for this simulation, please check %s.\n", fname);

			return -1;

		}

    	if (y1 < K.Ly / (- 2.0) || y1 > K.Ly / 2.0) 
		{
				printf("\nThe particles y positions are not in the expexpected range for this simulation, please check %s.\n", fname);

			return -1;

		}

		if (z1 < K.Lz / (- 2.0) || z1 > K.Lz / 2.0) 
		{
				printf("\nThe particles z positions are not in the expexpected range for this simulation, please check %s.\n", fname);

			return -1;

		} //TODO find a more readable way 
 
	    molecule[i].x = side_x_minus1 * x1;
        molecule[i].y = side_y_minus1 * y1; 
        molecule[i].z = side_z_minus1 * z1; 

        count_scan += fread(&molecule[i].e_x, sizeof(double), 1, fp);
        count_scan += fread(&molecule[i].e_y, sizeof(double), 1, fp);
        count_scan += fread(&molecule[i].e_z, sizeof(double), 1, fp);
    } 


    fclose(fp);

	if(count_scan != (K.PartNum * kDIM * 2) )
	{
		printf("\nFile does not contain 6 parameters (3 positions and 3 directions), please check %s.\n", fname);
		printf("\nParameters contained =  %u \n", count_scan);

			return -1;
	}
	fprintf (stderr, "\n%s read.\n", fname);
	
	sprintf(fname,"%s/turbField-%d.pos",work_dir, next_loop);
    fp = fopen(fname, "r");
 	if ( !(fp = fopen(fname, "r") ) ) 
		{
			fprintf(stderr, "\nError opening file %s;  Program aborted!\n\n", fname);
			return -1;

		}

   count_scan = 0;

    for (unsigned i = 0; i < K.PartNum; ++i)
    {
        count_scan += fread(&turb_field[i].e_x, sizeof(double), 1, fp);
        count_scan += fread(&turb_field[i].e_y, sizeof(double), 1, fp);
        count_scan += fread(&turb_field[i].e_z, sizeof(double), 1, fp);
    } 

//only used to reload    count_scan += fread(&last_t, sizeof(int), 1, fp);

    fclose(fp);
	if(count_scan != (K.PartNum * kDIM) )
	{
		fprintf (stderr, "\nFile does not contain 3 parameters, please check %s.\n", fname);
		fprintf (stderr, "\nParameters contained =  %u \n", count_scan);

		return -1;
	}
	fprintf (stderr, "\n%s read.\n", fname);
	return 0;
}



