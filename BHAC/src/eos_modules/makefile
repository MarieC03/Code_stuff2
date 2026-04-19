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
# Here we have the geometry related sources
#########################################################################
FOBJECTS += mod_eos$F mod_eos_polytrope$F mod_eos_tabulated$F mod_eos_tabulated_parameters$F mod_eos_idealgas$F mod_eos_hybrid$F
FOBJECTS += mod_eos_readtable_scollapse$F mod_eos_readtable_compose$F mod_eos_interpolation$F # mod_eos_inversion_T_to_eps$F
OBJECTS += mod_eos.o mod_eos_polytrope.o mod_eos_tabulated.o mod_eos_tabulated_parameters.o mod_eos_idealgas.o mod_eos_hybrid.o
OBJECTS += mod_eos_readtable_scollapse.o mod_eos_readtable_compose.o mod_eos_interpolation.o # mod_eos_inversion_T_to_eps.o

#########################################################################
# Dependencies:
#########################################################################
mod_eos.o : mod_eos_tabulated_parameters.o
mod_eos_tabulated.o : mod_eos_tabulated_parameters.o mod_eos_interpolation.o mod_eos.o mod_eos_readtable_scollapse.o mod_eos_readtable_compose.o mod_rootfinding.o
mod_eos_polytrope.o : mod_eos.o
mod_eos_idealgas.o : mod_eos.o
mod_eos_hybrid.o: mod_eos.o mod_eos_tabulated.o
mod_eos_readtable_scollapse.o : mod_eos.o mod_eos_tabulated_parameters.o # mod_eos_inversion_T_to_eps.o
mod_eos_readtable_compose.o : mod_eos.o mod_eos_tabulated_parameters.o
mod_eos_inversion_T_to_eps.o : mod_eos_tabulated_parameters.o
