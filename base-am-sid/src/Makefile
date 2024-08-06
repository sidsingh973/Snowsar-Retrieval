# Syntax Notes
#  "$^" refers to the output variable name
#  "$@" refers to the filenames of all prerequisites, separated by spaces
#  "-fcheck=all" does run-time checks on array bounds, args, etc.

#F90=gfortran-mp-4.6
F90=gfortran-9
CFLAGS=
LFLAGS=
MetroMEMLS : *.f90
	$(F90) $(CFLAGS) -o $@ $^

#clean:
#rm CommonVars.o MetroMEMLS.o CalcPex.o Invert_Matrix.o MCMC.o ObsModel.o RandNormal.o RandUniform.o 
#rm ReadHyperPar.o ReadObs.o ReadRunParams.o SetupMEMLS.o SS_Model.o Show_Theta.o
