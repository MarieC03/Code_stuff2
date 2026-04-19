#================================================================================
#
#    BHAC (The Black Hole Accretion Code) solves the equations of
#    general relativistic magnetohydrodynamics and other hyperbolic systems
#    in curved spacetimes.
#
#    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
#    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
#
#    This file is part of BHAC.
#
#    BHAC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BHAC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

#########################################################################
# Everything IO-related goes here:
#########################################################################
FOBJECTS += amrio$F  convert$F  kracken$F  analysis$F collapse$F slice$F \
            mod_shell$F mod_healpix_det$F
OBJECTS += amrio.o convert.o kracken.o  analysis.o collapse.o slice.o \
           mod_shell.o mod_healpix_det.o

#########################################################################
# Dependencies:
#########################################################################
amrio.o : kracken.o mod_forest.o slice.o collapse.o mod_shell.o mod_eos_hybrid.o mod_healpix_det.o
analysis.o : mod_constants.o
mod_shell.o : mod_interpolate.o convert.o
mod_healpix_det.o : mod_c2b_cart3d.o mod_post_process.o

