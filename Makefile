CC=gcc
FLAGS= -Wvector-operation-performance -Wshadow -Wsuggest-attribute=const -Wall -Wextra -march=native -ftree-loop-im -fopenmp -std=gnu99 #-fprofile-arcs -ftest-coverage #-funroll-loops 
ENDFLAGS= -lm -lgomp 
OBJS=main.o MDLoad.o MDConstants.o MDPostProcessing.o Turbulence.o
DEPS= MDLoad.h MDConstants.h Turbulence.h Molecule.h MDPostProcessing.h debug.h

EXE=post_proc
TESTDIR=tests/
ARGS= ../../testrun

ifeq ($(debug), 1)
RUNFLAGS= -O0 -g -Werror $(FLAGS)
else
RUNFLAGS= -O3 $(FLAGS)
endif
debug=0

.PHONY: all
all: $(EXE)

#debug and profiling
.PHONY: check
check: clean_tests $(EXE)
	cram -v $(TESTDIR)/*.t 

.PHONY: debug
debug: clean $(EXE)
	gdb --args ./$(EXE) $(ARGS) 

.PHONY: mem_check
mem_check: clean $(EXE)
	valgrind -v --leak-check=full --show-leak-kinds=all ./$(EXE) $(ARGS)

.PHONY: prof
prof: clean_prof $(EXE)
	valgrind --tool=callgrind ./$(EXE) $(ARGS)

.PHONY: time
time: clean $(EXE)
	time ./$(EXE) $(ARGS)	

#CLEAN
.PHONY: clean
clean: clean_prof clean_tests
	rm -f *.o $(EXE) 

.PHONY: clean_tests
clean_tests:
	rm -f $(TESTDIR)/*.err

.PHONY: clean_prof
clean_prof:
	rm -f gmon.out prof *.gcov *gcno

.PHONY: clean_sim
clean_sim:
	rm -f *.S output.dat *.pos

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



