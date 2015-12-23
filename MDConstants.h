#ifndef MDCONSTANTS_H
#define MDCONSTANTS_H

#include <stdio.h> // FILE*

#define kDIM  3

typedef struct MDConstants
{
    //input.dat consts
    //indentation represents actual file layout
    unsigned iteration_num;
    unsigned N[kDIM];
    double L[kDIM];
    double dcut;
    char starting_branch[50];
    double v_0; 
    double gamma_tilde; 
    double delta_t; 
    double gamma_rel; 
    double kappa;   unsigned NF;
    unsigned t_gap, deltaS;    

    //derived consts
    unsigned PartNum;
    unsigned SnapshotNum; 
    double side_minus1[kDIM];
//TODO better name for side_minus1

} MDConstants;

int
Initialize (MDConstants* MDConsts, char argv[]);

extern char work_dir[];

//helper functions
void
PrintVals
(
    const char message[],
    double * v_ptr,
    unsigned dim
);

void
CountReadUnsigned
(
    unsigned* count_scan,
    FILE* fp,
    unsigned* v_ptr,
    const unsigned dim
);

void
CountReadDouble
(
    unsigned* count_scan,
    FILE* fp,
    double* v_ptr,
    unsigned dim
);


#endif /* MDCONSTANTS_H */
