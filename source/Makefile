#############################################################################################
# openLBMflow v1.0.1 Copyright (C) 2013 LBMflow
# Open Source Lattice Boltzmann Solver
# www.lbmflow.com
# open@lbmflow.com
#
# LICENSE
# The openLBMflow code is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# DISCLAIMER OF WARRANTY
# The code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# This software may contain errors that could cause failures or loss of data,
# and may be incomplete or contain inaccuracies.  You expressly acknowledge and agree
# that use of the openLBMflow software is at your sole risk.
# The openLBMflow software is provided 'AS IS' and without warranty of any kind.
#############################################################################################

CC     = gcc
CFLAGS = -O2 
LIBS   = -lm
DEPS   = 
OBJ    = openLBMflow.o 

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

openLBMflow:$(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o 

#############################################################################################
