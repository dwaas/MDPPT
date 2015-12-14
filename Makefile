CC=gcc
FLAGS= -Wvector-operation-performance -Wshadow -Wsuggest-attribute=const -Wall -Wextra -Werror -march=native -ftree-loop-im -fopenmp -pg -std=gnu99 #-fprofile-arcs -ftest-coverage #-funroll-loops 
ENDFLAGS= -lm -lgomp 
OBJS=main.o MDLoad.o MDConstants.o 
DEPS= MDLoad.h MDConstants.h Turbulence.h Molecule.h
EXE=post_proc
ARGS= ../../testrun

ifeq ($(debug), 0)
RUNFLAGS= -O3 $(FLAGS)
else
RUNFLAGS= -O0 -g $(FLAGS)
endif
debug=0

.PHONY: all
all: $(EXE)

.PHONY: debug
debug: clean new $(EXE)
	gdb --args ./$(EXE) $(ARGS) 

.PHONY: clean
clean: 
	rm -f *.o $(EXE) gmon.out prof *.gcov *gcno

.PHONY: new
new:
	rm -f *.S output.dat *.pos

#PROFILING and DEBUGGING
.PHONY: prof
prof:
	gprof $(EXE) gmon.out > prof


$(EXE): $(OBJS)
	$(CC) -o $@ $^ $(RUNFLAGS) $(ENDFLAGS)

%.o:  %.c $(DEPS) 
	$(CC) -c -o $@ $< $(RUNFLAGS)

#FLAGS=  -O3 -ipo -static  -c 
#MAINFLAG=  -O3 -ipo -static

#FLAGS= -openmp -parallel -O3 -ipo -static -c 
#MAINFLAG= -openmp -parallel -O3 -ipo -static

#FLAGS= -openmp -parallel -fast  -c
#MAINFLAG= -openmp -parallel -fast 



