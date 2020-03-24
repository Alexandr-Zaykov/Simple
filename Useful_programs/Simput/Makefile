#
# ifort compilation of program Simput)
#
FORTS = $(wildcard *.f90)
OBJS  = $(FORTS:%.f90=%.o)
ifndef ($(prefix))
        prefix = ~/bin
endif
Target= $(prefix)/Simput
.SUFFIXES: .f90


all: $(Target)

$(Target): $(OBJS)
	ifort -openmp -O3 -mkl -o $(Target) $(OBJS)

.f90.o:
	ifort -openmp -O3 -c $<

clean:
	rm -f *.o

veryclen:
	rm -f *.o $(Target)

