ALL: heat_int1d


heat_int1d: u_icbc.o create_bintree.o chebexps.o prini.o cfgt.o gauss1d.o dpotentials.o heat_int1d.o int1d_dr.o lu.o


	gfortran -o heat_int1d $(FFLAGS) u_icbc.o create_bintree.o chebexps.o prini.o cfgt.o gauss1d.o dpotentials.o heat_int1d.o int1d_dr.o lu.o


.f.o:
	gfortran -c $(FFLAGS) $*.f
