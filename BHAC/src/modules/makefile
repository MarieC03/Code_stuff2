#########################################################################
# Here we have some auxilary modules
#########################################################################
FOBJECTS += mod_constants$F   mod_timing$F F              \
mod_odeint$F  mod_oneblock$F mod_lu$F mod_rtsafe2D$F	  \
mod_interpolate$F mod_ngmres$F mod_random$F               \
mod_rootfinding$F mod_lnsrch$F

OBJECTS += mod_constants.o    mod_timing.o                \
mod_odeint.o  mod_oneblock.o mod_lu.o mod_rtsafe2D.o	  \
mod_interpolate.o mod_ngmres.o mod_random.o               \
mod_rootfinding.o mod_lnsrch.o


#########################################################################
# Dependencies:
#########################################################################
mod_metric.o : mod_lu.o
mod_ngmres.o: mod_lu.o
mod_rootfinding.o: mod_lu.o mod_lnsrch.o
mod_interpolate.o: slice.o mod_forest.o
