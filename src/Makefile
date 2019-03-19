NC_LINK_CMD := $(shell nf-config --flibs)
NC_LINK_CMD += $(shell nc-config --libs)
NC_INC_CMD  := -I$(shell nc-config --includedir) -I$(shell nf-config --includedir)

LINK :=$(NC_LINK_CMD)
INCLUDE :=$(NC_INC_CMD)

ifeq ($(FC),mpif90)
  TEMP_VAR := $(shell mpif90 -show ) 
  FC := $(word 1,$(TEMP_VAR))
endif
ifeq ($(FC),)
  FC := ifort
endif
ifeq ($(FC),ifort)
  FCFLAGS := -cpp -w -auto -noalign -fautomatic -DMDL_BOX -DUSE_REAL8 -traceback
  # Debug
  ifdef DEBUG
    FCFLAGS += -g -check bounds -O0
    FCFLAGS += -fpe0 -ftrapuv -check arg_temp_created
    FCFLAGS += -finit-real=snan
    FLFLAGS := -g
  else
    # Performance
    FCFLAGS += -O2 -vec-report0
    FLFLAGS :=
  endif
  FREEFORM := -free
  R8       := -r8
endif
ifeq ($(FC),gfortran)
  FCFLAGS := -cpp -w -std=legacy -fautomatic -DMDL_BOX -DUSE_REAL8 -fbacktrace
  # Debug
  ifdef DEBUG
    FCFLAGS += -g -fbounds-check
    FCFLAGS += -ffpe-trap=invalid,zero,overflow
    FCFLAGS += -finit-real=snan
    FLFLAGS := -g
  else
    # Performance
    FCFLAGS += -O3 -funroll-loops
    FLFLAGS :=
  endif
  FREEFORM := -ffree-form -ffree-line-length-none
  R8       := -fdefault-real-8 -fdefault-double-8
endif

PROGRAM := sect_aer
PROJ_OBJ := $(PROGRAM).o

LD   :=$(FC) $(FCFLAGS) $(INCLUDE)

SRC  :=$(wildcard *.F) $(wildcard *.F90)
TMP  :=$(SRC:.F=.o)
OBJS :=$(TMP:.F90=.o)

MODS :=$(OBJS:.o=.mod)

default: exe

# Compiler steps for all objects
%.o %.mod : %.F
	$(LD) $(LINK) -c $<
%.o %.mod : %.F90
	$(LD) $(LINK) -c $(FREEFORM) $<

# Linker
exe: $(OBJS)
	$(LD) $(LINK) $(OBJS) $(FLFLAGS) -o $(PROGRAM)
	cp $(PROGRAM) ../bin/.

debug:
	@echo "SRC  = $(SRC)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "PROGRAM = $(PROGRAM)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(MODS)
	rm -rf ../bin/*

.PHONY: debug default clean

# Dependencies
physconstants.o                   : physconstants.F precision_mod.o

precision_mod.o                   : precision_mod.F

sect_aer_data_mod.o               : sect_aer_data_mod.F90 physconstants.o precision_mod.o

sect_aer_mod.o                    : sect_aer_mod.F90 sect_aer_data_mod.o physconstants.o precision_mod.o sect_aux_mod.o

sect_aux_mod.o                    : sect_aux_mod.F90 precision_mod.o

main_sect.o                       : main_sect.F90 sect_aer_mod.o sect_aer_data_mod.o precision_mod.o physconstants.o sect_aux_mod.o
