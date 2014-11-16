FF = gfortran
INCLUDE = -I/opt/local/include
LINK = -L/Users/Adem/Dropbox/lapack-3.5.0/
FLAGS = -llapack -lrefblas
EXTRAFLAGS = -g -ffree-line-length-0 -fopenmp
objects = array_work_asymm.o geometry_and_data_asymm.o nest_multy_mode_asymm.o beam_asymm.o field_asymm.o rab.o fft_mod.o

execute_assym : $(objects)
	$(FF) -o execute_assym $(objects) $(LINK) $(FLAGS) $(EXTRAFLAGS)

array_work_asymm.o : array_work_asymm.f90
	$(FF) -c array_work_asymm.f90 $(EXTRAFLAGS)

geometry_and_data_asymm.o : array_work_asymm.o rab.o geometry_and_data_asymm.f90
	$(FF) -c geometry_and_data_asymm.f90 $(EXTRAFLAGS)

nest_multy_mode_asymm.o : array_work_asymm.o geometry_and_data_asymm.o nest_multy_mode_asymm.f90 beam_asymm.o field_asymm.o rab.o fft_mod.o
	$(FF) -c nest_multy_mode_asymm.f90 $(EXTRAFLAGS)

beam_asymm.o : array_work_asymm.o rab.o beam_asymm.f90
	$(FF) -c beam_asymm.f90 $(EXTRAFLAGS)

field_asymm.o : array_work_asymm.o rab.o field_asymm.f90
	$(FF) -c field_asymm.f90 $(EXTRAFLAGS)

rab.o : array_work_asymm.o rab.f90
	$(FF) -c rab.f90 $(EXTRAFLAGS)

fft_mod.o : fft_mod.f90
	$(FF) -c fft_mod.f90 $(EXTRAFLAGS)

clean :
	rm *.o *.mod execute_assym
