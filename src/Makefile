FC = gfortran
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
PROGRAM := sect_aer
PROJ_OBJ := $(PROGRAM).o

FREEFORM := -ffree-form -ffree-line-length-none
R8       := -fdefault-real-8 -fdefault-double-8

LINK :=
LD   :=$(FC) $(FCFLAGS)

SRC  :=$(wildcard *.F) $(wildcard *.F90)
TMP  :=$(SRC:.F=.o)
OBJS :=$(TMP:.F90=.o)

MODS :=$(OBJS:.o=.mod)

default: exe

# Compiler steps for all objects
%.o %.mod : %.F
	$(LD) -c $<
%.o %.mod : %.F90
	$(LD) -c $(FREEFORM) $<

# Linker
exe: $(OBJS)
	$(LD) $(OBJS) $(FLFLAGS) -o $(PROGRAM)
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

sect_aer_mod.o                    : sect_aer_mod.F90 sect_aer_data_mod.o physconstants.o precision_mod.o

main_sect.o                       : main_sect.F90 sect_aer_mod.o sect_aer_data_mod.o precision_mod.o physconstants.o sect_aux_mod.o

sect_aux_mod.o                    : sect_aux_mod.F90 precision_mod.o
