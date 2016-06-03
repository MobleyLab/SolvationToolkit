#!/usr/bin/env python
######################################################################
# SolvationToolkit: A toolkit for setting up molecular simulations of mixtures
# Copyright 2011-2015 UC Irvine and the Authors
#
# Authors: David Mobley and Gaetano Calabro
# With thanks to Kyle Beauchamp, whose liquid_tools.py provided an initial basis for this code
# (https://github.com/choderalab/LiquidBenchmark/blob/master/src/simulation/liquid_tools.py) in April 2015
#
#This library is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2.1 of the License, or (at your option) any later version.
#
#This library is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with this library; if not, see <http://www.gnu.org/licenses/>.
######################################################################


from solvationtoolkit.solvated_mixtures import *
from openeye.oechem import *
from openeye.oeiupac import *


liquids = MixtureSystem()

liquids.addComponent(mole_fraction=0.2, name='ethane')
liquids.addComponent(mole_fraction=0.0, name='toluene')
liquids.addComponent(mole_fraction=0.4, name='methane')
liquids.addComponent(name='water')

liquids.build(amber=True, gromacs=True)


