openLBMflow v1.0.1 Copyright (C) 2013 LBMflow
=============================================
Open Source 2D/3D single and multiphase Lattice Boltzmann Code
--------------------------------------------------------------

The most recent release of openLBMflow can be downloaded at http://www.lbmflow.com/ or to the Email open@lbmflow.

Version 1.0.1
Release Date: 06 August 2013

CONTENTS:
	1. What is openLBMflow
	2. Installation
	3. Compilation
	4. How to configure model
	5. License
	6. Disclaimer of warranty

1. WHAT IS openLBMflow
	The openLBMflow is an open source Lattice Boltzman solver.
	Capable to simulate 2D or 3D, single or multiphase fluid flow.
	Single phase is a standard BGK lattice boltzman model,
	and Multiphase is a popular Schan-Chen BGK model.
	2D version is implemented on D2Q9 lattice and 3D on D3Q19 lattice.

2. INSTALLATION
	Unpack zip file to a folder on your local hard drive.
	You can compile the source code or run the binary file from bin/<system> folder.
	The precompiled binary files are available for windows and linux.
	To compile examples files replace openLBMFlow_conf.c configuration file in
	the source directory with the openLBMFlow_conf.c from examples directory,
	and compile the code with new configuration file.
	The program will create an output folder with data saved in VTK image format,
	to view results open .pvd file in Paraview software. This will load all
	saved timesteps of your simulation.

3. COMPILATION
	Command Line:
	To compile the solver from command line use:
	gcc -O2 -o openLBMflow openLBMflow.c -lm

	Linux:
	Under Linux just use make function in shell to compile the code
	according to Makefile provided in this release.

	Windows:
	Under Windows, you can install mingw GNU C/C++ compiler
	and use make.bat file to compile executable file.
	Alternatively you can install any IDE (i.e. Code:Blocks, NetBins, CodeLitle, VC...).

	Warning
	If you have any problems with compilation process,
	make sure that you have installed an C compiler in your system.
	if you are using the intel compiler use cc instead of gcc.

4. HOW TO CONFIGURE MODEL
	The model configuration can be changed in openLBMFlow_conf.c file.
	Default configuration is 3D Schan-Chen Multiphase model with 2 droplets.
	The droplets coalesce and due to gravity pointing down the formed droplet
	is falling down on a solid surface.
	By tuning the 'rho_boundary' parameter, static contact angle
	and wettability of surface can be changed.
	Several examples are provided in order to show potential of the solver.
	All examples are in 3D, but you can turn them into 2D simulation
	by deleting '#define Lattice3D ' line in openLBMFlow_conf.c file.

5. LICENSE
	The openLBMflow code is a free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

6. DISCLAIMER OF WARRANTY
	The code is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	This software may contain errors that could cause failures or loss of data,
	and may be incomplete or contain inaccuracies.  You expressly acknowledge and agree
	that use of the openLBMflow software is at your sole risk.
	The openLBMflow software is provided 'AS IS' and without warranty of any kind.
