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
FOBJECTS += mod_transform$F geometry$F mod_metric$F mod_bl$F mod_coord_cfc$F
FOBJECTS += mod_metric_aux$F
OBJECTS += mod_transform.o geometry.o mod_metric.o mod_bl.o mod_coord_cfc.o
OBJECTS += mod_metric_aux.o

#########################################################################
# Dependencies:
#########################################################################
mod_metric.o : mod_metric_aux.o mod_transform.o mod_rtsafe2D.o mod_interpolate.o mod_cfc_parameters.o mod_coord_cfc.o
mod_coord_cfc.o: mod_metric_aux.o 
geometry.o initialize_vars.o : mod_metric.o mod_bl.o mod_limiter.o
mod_transform.o : mod_bl.o
