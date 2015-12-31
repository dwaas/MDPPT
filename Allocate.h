#ifndef ALLOCATE_H
#define ALLOCATE_H

#include <stdlib.h> // malloc, free

#include "debug.h"

#define FREE(p)										\
			free (p);								\
			p = NULL; 

#define FREE2D(p,size_x)							\
			for (unsigned i = 0; i < (size_x); ++i)	\
			{ 										\
				FREE (p[i]);						\
			}										\
			FREE (p);

#define FREE3D(p,size_y,size_x)						\
			for (unsigned j = 0; j < (size_x); ++j)	\
			{ 										\
				FREE2D (p[j], size_y);				\
			}										\
			FREE(p);								

#define ALLOC(p,size_x)								\
			if (p) free (p);						\
			p = calloc ((size_x), sizeof(*p));		\
			check_mem (p);

#define ALLOC2D(p,size_x,size_y)					\
			ALLOC(p,(size_x));						\
			for (unsigned i = 0; i < (size_x); ++i)	\
			{										\
				ALLOC (p[i], (size_y) );			\
			}																		

#define ALLOC3D(p,size_x,size_y,size_z)				\
			ALLOC (p, (size_x) );					\
			for (unsigned j = 0; j < (size_x); ++j)	\
			{										\
				ALLOC2D (p[j], (size_y), (size_z) );\
			}													

#endif /* ALLOCATE_H */
