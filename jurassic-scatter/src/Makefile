# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
#EXC = brightness climatology collect formod interpolate kernel limb nadir planck raytrace retrieval tab2bin
EXC = formod

# Library directories...
#LIBDIR = -L /usr/lib64 
LIBDIR = 
LIBDIR += -L /p/software/jusuf/stages/2020/modules/all/Compiler/GCC/9.3.0
LIBDIR += -L /p/software/jusuf/stages/2020/software/CUDA/11.0/lib64
LIBDIR += -L ../jurassic-gpu/src


# Include directories...
#INCDIR = -I /usr/local/include
INCDIR += -I ../common-header
INCDIR += -I ../jurassic-gpu/src/inter_folder

# OpenMP...
OMP = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler with/without openMP and without MPI...
# -----------------------------------------------------------------------------

# Compiler...
#CC = gcc

# CFLAGS...
#CFLAGS = $(INCDIR) -pedantic-errors -Wall -Wextra -Werror -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long -Winline -fno-common -DHAVE_INLINE

# Flags for OpenMP...
# ifdef OMP
# CFLAGS += -fopenmp
# endif

# -----------------------------------------------------------------------------
# Set flags for open MPI C wrapper compiler with/without openMP...
# -----------------------------------------------------------------------------
# # Compiler...
CC = mpicc 

# # CFLAGS...
CFLAGS = $(INCDIR) -pedantic-errors -Wall -Wextra -Werror -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long -Winline -fno-common -DHAVE_INLINE -DMPI

# # Flags for OpenMP...
ifdef OMP
CFLAGS += -fopenmp
endif
# -----------------------------------------------------------------------------

CFLAGS += -D WORK_QUEUE


# Flags for profiling...
ifdef PROF
CFLAGS += -O2 -g -pg -ggdb3
else
CFLAGS += -O3
endif

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm -lstdc++ -lcudart -lutil

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)
	rm -f *~

$(EXC): %: %.c  atmosphere.o control.o sca_forwardmodel.o workqueue.o sca_scatter.o Makefile
	$(CC) $(CFLAGS) -o $@ $< atmosphere.o control.o sca_forwardmodel.o  workqueue.o sca_scatter.o  $(LDFLAGS)

atmosphere.o: atmosphere.c 
	$(CC) $(CFLAGS) -c atmosphere.c

control.o: control.c
	$(CC) $(CFLAGS) -c control.c

sca_forwardmodel.o: sca_forwardmodel.c 
	$(CC) $(CFLAGS) -c sca_forwardmodel.c

#misc.o: misc.c 
#	$(CC) $(CFLAGS) -c misc.c

#retrievalmodel.o: retrievalmodel.c 
#	$(CC) $(CFLAGS) -c retrievalmodel.c

sca_scatter.o: sca_scatter.c 
	$(CC) $(CFLAGS) -c sca_scatter.c

workqueue.o: workqueue.c 
	$(CC) $(CFLAGS) -c workqueue.c



bak:
	tar czf ../bak/jurassic_`date +"%y%m%d%H%M"`.tgz Makefile *.c *.h

clean:
	rm -f $(EXC) *.o *~
