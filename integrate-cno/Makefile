# Make the HSE Integrator

SUNINCDIR = /home/eugene/local/sundials/instdir/include
SUNLIBDIR = /home/eugene/local/sundials/instdir/lib

F90 	= gfortran
FFLAGS 	= -O3 -ffree-line-length-none #-ffpe-trap=invalid,zero,overflow
SUNINCLUDE = -I ${SUNINCDIR}
INCLUDE = 
F90LINK	= gfortran
CC	= gfortran

LLIBSS	= -L${SUNLIBDIR} -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial
LLIBSM	= -L/usr/lib64 -llapack -lblas

SRCS 	= cvodenet.f90 cno_nuc9.f90 integrator.f90 physical_constants.f90 data_wrangler.f90 parameters.f90 cvode_parameters.f90 
OBJS	= cvodenet.o cno_nuc9.o cvode_parameters.o data_wrangler.o integrator.o parameters.o physical_constants.o 

MAIN	= intnet

all:	${MAIN}
	@echo	Compilation complete!

%.mod : %.o
	@echo Compiled module.

call_net : cno_nuc9.mod physical_constants.mod
	${F90} ${FFLAGS} call_net.f90 cno_nuc9.o physical_constants.o -o call_net

integrator.o : integrator.f90 cvodenet.mod cno_nuc9.mod data_wrangler.mod parameters.mod cvode_parameters.mod
	${F90} ${FFLAGS} -c $< -o $@

cvode_parameters.o : cvode_parameters.f90
	${F90} ${FFLAGS} -c $< -o $@

cno_nuc9.o : cno_nuc9.f90 physical_constants.mod
	${F90} ${FFLAGS} -c $< -o $@

cvodenet.o : cvodenet.f90 cno_nuc9.mod physical_constants.mod cvode_parameters.mod
	${F90} ${FFLAGS} -c $< -o $@

data_wrangler.o : data_wrangler.f90
	${F90} ${FFLAGS} -c $< -o $@

parameters.o : parameters.f90 cvode_parameters.mod cno_nuc9.mod
	${F90} ${FFLAGS} -c $< -o $@

physical_constants.o : physical_constants.f90
	${F90} ${FFLAGS} -c $< -o $@

$(MAIN): ${OBJS}
	@echo '${CC} ${INCLUDE} -o ${MAIN} ${OBJS} ${LLIBSS} ${LLIBSM}'
	${CC} ${INCLUDE} -o ${MAIN} ${OBJS} ${LLIBSS} ${LLIBSM}

clean:
	rm cvodenet.o  cvode_indices.o  cvode_parameters.o  data_wrangler.o  integrator.o  parameters.o  physical_constants.o call_net *.mod *~ ${MAIN}

depend: ${SRCS}
	makedepend ${INCLUDE} $^
