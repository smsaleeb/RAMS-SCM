#Makefile

###############################################################################
# Include definitions
###############################################################################
MAKE=/usr/bin/make
MODEL=./src

#F_COMP=/home/smsaleeb/intel/composer_xe_2011_sp1.8.273/bin/intel64/ifort

F_COMP=/home/smsaleeb/software/mpich-3.3.2/bin/mpif90
PAR_LIBS=-L/home/smsaleeb/software/mpich-3.3.2/lib -lmpifort -lmpi

F_OPTS=-free -O2 -fp-model precise
LOADER_OPTS= -free -O2 -fp-model precise

LIBS=-L/usr/lib/x86_64-linux-gnu -lrt
ARCH=ar rsU

###############################################################################
# Compiler commands
###############################################################################
F_CD = $(F_COMP) -c $(F_OPTS)
F_CM = $(F_CD) $< && $(ARCH) $@ $(<F:.f90=.o) && rm -f $(<F:.f90=.o)

################################################################################
## File extension rules
################################################################################

$(ARC)(%.o): %.f90; $(F_CM)

################################################################################
# Define objects
################################################################################
OBJ = $(ARC)($(MODEL)/lib/rconstants.o) \
      $(ARC)($(MODEL)/memory/grid_dims.o) \
      $(ARC)($(MODEL)/micro/micphys.o) \
      $(ARC)($(MODEL)/io/io_params.o) \
      $(ARC)($(MODEL)/memory/mem_other.o) \
      $(ARC)($(MODEL)/memory/mem_basic.o) \
      $(ARC)($(MODEL)/memory/mem_leaf.o) \
      $(ARC)($(MODEL)/memory/mem_micro.o) \
      $(ARC)($(MODEL)/memory/mem_radiate.o) \
      $(ARC)($(MODEL)/mpi/node_mod.o) \
      $(ARC)($(MODEL)/memory/mem_grid.o) \
      $(ARC)($(MODEL)/lib/numutils.o) \
      $(ARC)($(MODEL)/lib/dateutils.o) \
      $(ARC)($(MODEL)/lib/therm_lib.o) \
      $(ARC)($(MODEL)/bc/rbnd_trsetsns.o) \
      $(ARC)($(MODEL)/core/rthrm.o) \
      $(ARC)($(MODEL)/core/rtimi.o) \
      $(ARC)($(MODEL)/micro/aero_include.o) \
      $(ARC)($(MODEL)/radiate/rrad3.o) \
      $(ARC)($(MODEL)/radiate/rad_driv.o) \
      $(ARC)($(MODEL)/radiate/rad_aero.o) \
      $(ARC)($(MODEL)/radiate/rad_mclat.o) \
      $(ARC)($(MODEL)/radiate/rrad2.o) \
      $(ARC)($(MODEL)/micro/mic_aero.o) \
      $(ARC)($(MODEL)/micro/mic_chknan.o) \
      $(ARC)($(MODEL)/micro/mic_init.o) \
      $(ARC)($(MODEL)/micro/mic_init_scm.o) \
      $(ARC)($(MODEL)/micro/mic_coll.o) \
      $(ARC)($(MODEL)/micro/mic_driv.o) \
      $(ARC)($(MODEL)/micro/mic_misc.o) \
      $(ARC)($(MODEL)/micro/mic_adj.o) \
      $(ARC)($(MODEL)/micro/mic_vap.o) \
      $(ARC)($(MODEL)/micro/mic_nuc.o) \
      $(ARC)($(MODEL)/micro/mic_nuctab.o) \
      $(ARC)($(MODEL)/micro/mic_nucpre.o) \
      $(ARC)($(MODEL)/micro/mic_tabs.o) \
      $(ARC)($(MODEL)/micro/mic_gamma.o) \
      $(ARC)($(MODEL)/micro/aero_sources.o) \
      $(ARC)($(MODEL)/micro/aero_deposit.o) \
      $(ARC)($(MODEL)/micro/write_scm.o)
###############################################################################
# Define archive and executable names
###############################################################################
BASE=rams-scm
EXE=$(BASE)
ARC=$(BASE).a

# Define main source.
MAIN_OBJ = ./main.o
MAIN = $(MODEL)/main.f90

# Define targets.

all: $(EXE)

$(EXE): $(ARC) $(MAIN) FORCE
	@echo ""
	$(F_COMP) -o $(EXE) $(MAIN_OBJ) $(LOADER_OPTS) $(ARC) $(LIBS) $(PAR_LIBS)
	rm -f *.o
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAIN): FORCE
	@echo ""
	$(F_CD) $@

$(ARC): $(OBJ)

FORCE:

check: FORCE
	@echo ""
	check

clean:
	@echo ""
	rm -f $(ARC) $(EXE) $(BASE) *.o *.mod *.f
	@echo ""
