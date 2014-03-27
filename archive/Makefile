# -*- makefile -*-

F2PY = f2py
F2PY_FLAGS = --opt="-O3" --fcompiler=gfortran
F2PY_LIBS = -L/Users/brentminchew/bin/lapack-3.4.2 -llapack -lrefblas -lgfortran

TARGETS = ./_reconutils

all: $(TARGETS)

./_reconutils: ./_reconutils.f90
	 $(F2PY) $(F2PY_FLAGS) $(F2PY_LIBS) -c $< -m $@

clean:
	 rm -f $(TARGETS).so

again: clean all 

