#F90 = lf95
#FFLAGS =  -O --nap --tpp --ntrace
#FLAGS = --ap --chk -g  --warn --trap
#F90 = lf95
#F90FLAGS = --ap --chk -g  --warn --trap
#LD = $(F90)
#LDFLAGS
#F90 = gfortran 
F90 = mpif90
#F90 = ifort
#F90 = gfortran -ffixed-line-length-132
#FFLAGS = -g -w90 -w95 -align dcommons
#F90OPT = -O3 -align dcommons -tpp7
#F90OPT = -mcmodel=medium -shared-intel -O2  -C
#F90OPT = -g -C -traceback -warn 
#F90OPT = -O3 
F90OPT = -g -fbacktrace -fbounds-check
#F90OPT = -g -C #-Wall -fbacktrace
libdir       = /home/pierre/sundials-2.6.2/sundials/install/lib
LIBRARIES = -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial ${LIBS}
LINKFLAGS = -Wl,-rpath,/home/pierre/sundials-2.6.2/sundials/install/lib

.PRECIOUS:

.SUFFIXES: .f .f90 .o .a

all: Couplage 

clean:
	@rm -f *~ *.o *.mod *__genmod.f90 Couplage

Couplage: prec.o parametres.o tools.o fonctioncvode.o stochastique.o main.o 
	$(F90) $(F90OPT) -o Couplage prec.o parametres.o tools.o fonctioncvode.o stochastique.o main.o -L${libdir} ${LIBRARIES}  ${LINKFLAGS}
prec.o: prec.f90
	$(F90) $(F90OPT) -c prec.f90 -L${libdir} ${LIBRARIES} ${LINKFLAGS}
parametres.o: prec.f90 parametres.f90
	$(F90) $(F90OPT) -c parametres.f90 -L${libdir} ${LIBRARIES} ${LINKFLAGS}
fonctioncvode.o: prec.f90 parametres.f90 fonctioncvode.f90 
	$(F90) $(F90OPT) -c fonctioncvode.f90 -L${libdir} ${LIBRARIES} ${LINKFLAGS}
tools.o: prec.f90 parametres.f90 tools.f90 
	$(F90) $(F90OPT) -c tools.f90 -L${libdir} ${LIBRARIES} ${LINKFLAGS}
stochastique.o: prec.f90 parametres.f90 tools.f90 stochastique.f90 
	$(F90) $(F90OPT) -c stochastique.f90 -L${libdir} ${LIBRARIES} ${LINKFLAGS}
main.o: prec.f90 parametres.f90 tools.f90 fonctioncvode.f90 stochastique.f90 main.f90 
	$(F90) $(F90OPT) -c main.f90 -L${libdir} ${LIBRARIES} ${LINKFLAGS} 
