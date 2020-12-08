# use f77 to compile - most strict fortran compiler I could find on Ubunutu that was easy to install
# -freal-4-real-8 ensures that double precision is used exclusively. 
FF = f77
FFLAGS = -c -Wall -freal-4-real-8 -std=legacy

#Objects allows multiple experiments to be created at once
objects = LoopExperiments.out
all : $(objects)

#links the program to the base programs 
%.out : FDVM2-Solver.o InitialConditionProblems.o MiscellaneousFunctions.o %.a
	$(FF) $^ -o $@

#compiles the program that relies on and must be linked to base programs
%.a : %.f FDVM2-Solver.o InitialConditionProblems.o MiscellaneousFunctions.o
	$(FF) $(FFLAGS) $< -o $@

%.o : %.f
	$(FF) $(FFLAGS) $< -o $@

clean:
	rm -f ./*.o *~ ./*.out

fresh : clean all
