# Set the compiler flags
#FFLAGS = -fast#-fast #-arch pn4 -tpp7 -tune pn4 -ax
#FFLAGS = -r8 -fast -fastsse -O3 -Mipa=fast -Minline -Munroll	#PGI optimised
#FFLAGS = -Mbounds -g 				#PGI Debug
#FFLAGS = -O3 -fast                            	#Intel
#FFLAGS = -fpe0 -nothreads -traceback -fltconsistency -CB -g -inline_debug_info #Intel Debug

FFLAGS = -O1

# Set some of the build parameters
TARGET = lare2d

# Uncomment one of the following lines if on a cluster.
#COPSON = -compiler intel
#MHDCLUSTER = -f90=pgf90 -DMHDCLUSTER -fpic

#Uncomment the following line to use Qmono viscosity
#QMONO = -DQ_MONO


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODULEFLAG = -module
MACHINEFLAGS = $(COPSON) $(MHDCLUSTER)
OPFLAGS = $(QMONO)
FC = mpif90 $(MACHINEFLAGS) $(OPFLAGS)
PREPROFLAGS = $(NONMPIIO)

OBJFILES = shared_data.o mpi_routines.o openboundary.o mpiboundary.o boundary.o normalise.o conduct.o diagnostics.o setup.o lagran.o  \
 remap.o xremap.o yremap.o zremap.o initial_conditions.o\
 output_cartesian.o iocontrol.o output.o iocommon.o input.o inputfunctions.o\
 input_cartesian.o eos.o neutral.o control.o\
 welcome.o lare2d.o
FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(OBJDIR):$(SRCDIR)/core:$(SRCDIR)/io/

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES))

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ $(SRCDIR)/core/*~ $(SRCDIR)/io/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: datatidy
datatidy:
	@rm -rf Data/*

.PHONEY: visit
visit:
	@cd VisIT;xml2makefile -clobber cfd.xml;make

.PHONEY: visitclean
visitclean:
	@cd VisIT;make clean;rm -f .depend
# All the dependencies
shared_data.o:shared_data.F90
mpi_routines.o:mpi_routines.f90 shared_data.o 
normalise.o:normalise.f90 shared_data.o
setup.o:setup.F90 shared_data.o normalise.o iocommon.o iocontrol.o input.o input_cartesian.o
mpiboundary.o: mpiboundary.f90 shared_data.o
openboundary.o: openboundary.f90 shared_data.o
boundary.o:boundary.f90 shared_data.o mpiboundary.o openboundary.o
xremap.o:xremap.f90 shared_data.o boundary.o
yremap.o:yremap.f90 shared_data.o boundary.o
zremap.o:zremap.f90 shared_data.o boundary.o
diagnostics.o:diagnostics.F90 shared_data.o boundary.o normalise.o output_cartesian.o output.o iocontrol.o eos.o
iocommon.o:iocommon.f90 shared_data.o
output.o:output.f90 shared_data.o iocommon.o
output_cartesian.o: output_cartesian.f90 shared_data.o iocommon.o output.o
iocontrol.o: iocontrol.f90 shared_data.o iocommon.o output.o input.o
input.o: input.f90 shared_data.o iocommon.o inputfunctions.o
inputfunctions.o: inputfunctions.f90 shared_data.o iocommon.o
conduct.o:conduct.f90 shared_data.o boundary.o eos.o
lagran.o:lagran.F90 shared_data.o boundary.o diagnostics.o conduct.o normalise.o eos.o neutral.o
remap.o:remap.f90 shared_data.o xremap.o yremap.o zremap.o
initial_conditions.o:initial_conditions.f90 shared_data.o normalise.o eos.o
eos.o:eos.F90 shared_data.o normalise.o
neutral.o: neutral.f90 shared_data.o boundary.o normalise.o eos.o
control.o: control.f90 shared_data.o normalise.o
welcome.o: welcome.f90 shared_data.o
lare2d.o:lare2d.f90 shared_data.o setup.o boundary.o diagnostics.o lagran.o remap.o mpi_routines.o welcome.o initial_conditions.o openboundary.o eos.o control.o