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
    unsigned Length;
    double InvLength;
    unsigned SnapshotNum; 

} MDConstants;


void Initialize (MDConstants* MDConsts, char argv[]);

#endif /* MDINPUTFILE_H */
