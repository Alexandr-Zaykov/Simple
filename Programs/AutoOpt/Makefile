#
# ifort compilation of program AutoOpt
#
FORTS = $(wildcard *.f90)
OBJS  = $(FORTS:%.f90=%.o)
ifndef ($(prefix))
        prefix = ~/bin
endif
Target= $(prefix)/AutoOpt
.SUFFIXES: .f90


all: $(Target)

$(Target): $(OBJS)
	ifort -O3 -mkl -openmp -o $(Target) $(OBJS)

.f90.o:
	ifort -O3 -openmp -c $<

clean:
	rm *.o

veryclen:
	rm *.o $(Target)

