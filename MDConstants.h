#ifndef MDINPUTFILE_H
#define MDINPUTFILE_H

typedef struct MDConstants
{
    //input.dat consts
    //indentation represents actual file layout
    unsigned iteration_num;
    unsigned Nx, Ny, Nz;
    double Lx, Ly, Lz;
    double dcut;
    char starting_branch[50];
    double v_0; 
    double gamma_tilde; 
    double delta_t; 
    double gamma_rel; 
    double kappa;   unsigned NF;
    unsigned t_gap, deltaS;    

    //derived consts
    unsigned DimNum;
    unsigned PartNum;
    double Length;
    unsigned SnapshotNum; 
    double InvLength;
    double 
            side_x_minus1, 
            side_y_minus1, 
            side_z_minus1;


} MDConstants;

void
Initialize (MDConstants* MDConsts, char argv[]);

extern char work_dir[];

#endif /* MDINPUTFILE_H */
