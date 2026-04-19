#########################################################################
# Here we have some auxilary modules
#########################################################################
FOBJECTS += handle_particles$F integrate_particles$F init_particles$F \
	    mod_particles$F mod_gridvars$F

OBJECTS += handle_particles.o integrate_particles.o init_particles.o  \
	    mod_particles.o mod_gridvars.o

#########################################################################
# Dependencies:
#########################################################################
mod_particles.o : mod_constants.o mod_gridvars.o slice.o
advance.o : mod_particles.o
amrvac.o : mod_gridvars.o
init_particles.o : mod_particles.o
