#
# ifort compilation of program SimOpt
#
FORTS = $(wildcard *.f90)
OBJS  = $(FORTS:%.f90=%.o)
USEMOD = $(shell grep -l 'use Declaration' *.f90)
ifndef ($(prefix))
	prefix = .
endif
Target= $(prefix)/SimOpt
.SUFFIXES: .f90 .mod

all: $(Target)

$(Target): $(OBJS)
	ifort -O3 -mkl -o $(Target) $(OBJS)

Declaration.o: Declaration.f90
	@touch $(USEMOD)
	ifort -O3 -c Declaration.f90

.f90.o:
	ifort -O3 -c $<

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod $(Target)

distclean: veryclean

$(filter-out Declaration.f90,$(FORTS)): declaration.mod
declaration.mod: Declaration.f90 Declaration.o
