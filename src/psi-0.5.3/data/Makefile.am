## @start 1
## Makefile.am -- Process this file with automake to produce Makefile.in
## @end 1
## Copyright (C) 2000 Gary V. Vaughan
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
## @start 1

EXTRA_DIST              = $(ENERGYDATA)
MAINTAINERCLEANFILES    = Makefile.in

ENERGYDATA              = \
pdb1afp.T.prd   pdb1afp.X.prd   pdb1afp._.p     pdb1g90.A.prd   pdb1v63.A.prd   xaa \
Makefile.in     pdb1afp.U.prd   pdb1afp.Y.prd   pdb1afp._.prd   pdb1hqb.A.prd   pdbtest.1.dat \
README          pdb1afp.V.prd   pdb1afp._.b     pdb1afp._.q     pdb1r48.A.dat   pdbtest.2.dat \
energy.int      pdb1afp.W.prd   pdb1afp._.o     pdb1g5m.A.prd   pdb1r48.A.prd   pdbtest.2.prd


##psidatadir               = $(srcdir)/data
psidatadir               = $(srcdir)

bin_SCRIPTS = data

data:
	@echo Making $@ in $(psidatadir)
	@cp $(psidatadir)/*.int $(top_builddir)/data 
	@cp $(psidatadir)/*.dat $(top_builddir)/data 
	@cp $(psidatadir)/*.prd $(top_builddir)/data 
	@cp $(psidatadir)/xaa $(top_builddir)/data 
	@cp $(psidatadir)/pdb1afp.* $(top_builddir)/data 

