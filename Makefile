# Specify a particular compiler with "make COMPILER=pgi", etc.
# Specify debugging flags with "make MODE=debug"
# Alternatively, these options can be specified as environment variables
# eg. "export COMPILER=gfortran" can be added to $HOME/.bashrc


# Compiler specific flags

# Note: you MUST specify a COMPILER option. None are specified by default.

ifeq ($(strip $(COMPILER)),)
  MAKECMDGOALS = error
error:
	@echo '*** ERROR ***'
	@echo 'You MUST set a value for the COMPILER variable'
	@echo ' eg. "make COMPILER=intel"'
	@echo 'Alternatively, you can add "export COMPILER=intel" to $$HOME/.bashrc'
endif

D = -D

# PGI
# ===
ifeq ($(strip $(COMPILER)),pgi)
  FFLAGS = -r8 -fast -fastsse -O3 -Mipa=fast,inline -Minfo # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -Mbounds -g                                     # Debug
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# Intel
# =====
ifeq ($(strip $(COMPILER)),intel)
  FFLAGS = -O3
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xHost # Optimised (B)
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xAVX  # Optimised (W)
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
             -C -g -heap-arrays 64 -warn -fpic
    ifeq ($(strip $(SYSTEM)),Darwin)
      FFLAGS += -Wl,-no_pie
    endif
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# gfortran
# ========
ifeq ($(strip $(COMPILER)),gfortran)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -Wall -Wextra -pedantic -fbounds-check \
             -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter \
						 -Wrealloc-lhs-all
    #FFLAGS += -ffpe-trap=underflow,denormal

    GNUVER := $(shell gfortran -dumpversion | head -1 \
        | sed 's/[^0-9\.]*\([0-9\.]\+\).*/\1/')
    GNUMAJOR := $(shell echo $(GNUVER) | cut -f1 -d\.)
    GNUMINOR := $(shell echo $(GNUVER) | cut -f2 -d\.)

    # gfortran-4.3
    GNUGE43 := $(shell expr $(GNUMAJOR) \>= 4 \& $(GNUMINOR) \>= 3)
    ifeq "$(GNUGE43)" "1"
      FFLAGS += -fbacktrace -fdump-core

      # gfortran-4.6
      GNUGE46 := $(shell expr $(GNUMINOR) \>= 6)
      ifeq "$(GNUGE46)" "1"
        FFLAGS += -Wno-unused-dummy-argument

        # gfortran-4.8
        GNUGE48 := $(shell expr $(GNUMINOR) \>= 8)
        ifeq "$(GNUGE48)" "1"
          FFLAGS += -Wno-target-lifetime
        endif
      endif
    endif
  endif
  MODULEFLAG = -I/usr/local/include -I$(OBJDIR) -J$(OBJDIR)
endif

# g95
# ========
ifeq ($(strip $(COMPILER)),g95)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g                                        # Debug
  endif
  MODULEFLAG = -fmod=$(OBJDIR)
endif

# IBM Bluegene
# ============
ifeq ($(strip $(COMPILER)),ibm)
  FFLAGS = -O5 -qhot -qipa # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -C -g -qfullpath -qinfo #-qkeepparm -qflttrap \
          -qnosmp -qxflag=dvz -Q! -qnounwind -qnounroll # Debug
    #FFLAGS = -O0 -qarch=qp -qtune=qp
    #FFLAGS = -qthreaded -qsmp=noauto -qsmp=omp # Hybrid stuff
  endif
  MODULEFLAG = -I$(OBJDIR) -qmoddir=$(OBJDIR)
  MPIF90 ?= mpixlf90_r

  # IBM compiler needs a -WF to recognise preprocessor directives
  D = -WF,-D
endif

# HECToR
# ========
ifeq ($(strip $(COMPILER)),hector)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -ea -ec -eC -eD -eI -en -hfp_trap -Ktrap=fp -m0 -M1438,7413
  endif
  MODULEFLAG = -em -I/usr/include -I$(OBJDIR) -J$(OBJDIR)
  MPIF90 ?= ftn
endif


MPIF90 ?= mpif90
FFLAGS += -I$(SDF)/include
FFLAGS += $(MODULEFLAG)
LDFLAGS = $(FFLAGS) -L$(SDF)/lib -lsdf

# Set some of the build parameters
TARGET = lare2d

# Set pre-processor defines
DEFINES := $(DEFINE)

# The following are a list of pre-processor defines which can be activated by 
# uncommenting some of the lines below.

# Uncomment to use Cauchy solution for predictor step B-field, othwerwise advective prediction
#DEFINES += $(D)CAUCHY

# Uncomment the following line to run with limiters on shock viscosity
#DEFINES += $(D)SHOCKLIMITER

# Uncomment the following line to allow shock viscosity in expanding cells
#DEFINES += $(D)SHOCKEXPANSION

# Uncomment the following line to run in single precision
#DEFINES += $(D)SINGLE

# Uncomment the following line to use first order scheme for resistive update
#DEFINES += $(D)FOURTHORDER

# Uncomment to add the 'nfs' file prefix required by some filesystems
#DEFINES += $(D)FILEPREFIX

# Don't generate any output at all. Useful for benchmarking.
#DEFINES += $(D)NO_IO

# Uncomment to enable the MPI error handler which is useful for debugging
#DEFINES += $(D)MPI_DEBUG


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

all: main

SDF := SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
  $(D)_MACHINE='"$(MACHINE)"'

SRCFILES = boundary.f90 conduct.f90 radiative.f90 control.f90 diagnostics.F90 \
  initial_conditions.f90 lagran.F90 lare2d.f90 mpi_routines.F90 \
  mpiboundary.f90 neutral.f90 normalise.f90 openboundary.f90 \
  random_generator.f90 remap.f90 setup.F90 shared_data.F90 version_data.F90 \
  welcome.f90 xremap.f90 yremap.f90 zremap.f90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(SDF)/src:$(OBJDIR)

-include $(SRCDIR)/COMMIT

ifeq ($(DONE_COMMIT),)
main: commit
else
main: $(FULLTARGET)
endif

# Rule to build the fortran files

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(FC) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

$(SDFMOD):
	$(MAKE) -C $(SDF)

clean:
	@rm -rf $(BINDIR) $(OBJDIR)

cleanall: tidy

tidy:
	@rm -rf $(OBJDIR) *~ *.pbs.* *.sh.* $(SRCDIR)/*~ *.log
	$(MAKE) -C $(SDF) cleanall

datatidy:
	@rm -rf Data/*

tarball:
	@sh $(SRCDIR)/make_tarball.sh

visit:
	@cd $(SDF)/../VisIt; ./build

visitclean:
	@cd $(SDF)/../VisIt; make clean; ./build -c; \
	  rm -rf .depend *.d *Info.C *Info.h CMake* cmake* Makefile

sdfutils:
	@cd $(SDF)/../C; make
	@cd $(SDF)/../utilities; sh build.sh

sdfutilsclean:
	@cd $(SDF)/../C; make clean
	@cd $(SDF)/../utilities; sh build.sh -c

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE

# All the dependencies

boundary.o: boundary.f90 mpiboundary.o random_generator.o shared_data.o
conduct.o: conduct.f90 neutral.o boundary.o shared_data.o
radiative.o: radiative.f90 boundary.o shared_data.o
control.o: control.f90 normalise.o shared_data.o
diagnostics.o: diagnostics.F90 boundary.o conduct.o shared_data.o \
  version_data.o $(SDFMOD)
initial_conditions.o: initial_conditions.f90 neutral.o diagnostics.o shared_data.o \
  boundary.o
lagran.o: lagran.F90 boundary.o conduct.o radiative.o neutral.o shared_data.o \
  openboundary.o remap.o
lare2d.o: lare2d.f90 boundary.o control.o diagnostics.o initial_conditions.o \
  lagran.o mpi_routines.o neutral.o normalise.o remap.o setup.o \
  shared_data.o welcome.o
mpi_routines.o: mpi_routines.F90 shared_data.o
mpiboundary.o: mpiboundary.f90 shared_data.o
neutral.o: neutral.f90 boundary.o shared_data.o
normalise.o: normalise.f90 shared_data.o
openboundary.o: openboundary.f90 shared_data.o
random_generator.o: random_generator.f90
remap.o: remap.f90 boundary.o shared_data.o xremap.o yremap.o zremap.o
setup.o: setup.F90 diagnostics.o shared_data.o version_data.o welcome.o \
  $(SDFMOD)
shared_data.o: shared_data.F90 $(SDFMOD)
version_data.o: version_data.F90 COMMIT
welcome.o: welcome.f90 shared_data.o version_data.o
xremap.o: xremap.f90 boundary.o
yremap.o: yremap.f90 boundary.o
zremap.o: zremap.f90 boundary.o
