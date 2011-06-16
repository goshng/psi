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

EXTRA_DIST              = $(TESTS)
MAINTAINERCLEANFILES    = Makefile.in

testsubdir              = testSubDir

TESTS_ENVIRONMENT       = top_builddir=$(top_builddir)

TESTS                   = psi seq grid energy rng gibbs mcmc ds bf mc nsv \
                          sim

distclean-local:
	-rm -rf $(testsubdir)
## @end 1

INCLUDES	= -I$(top_builddir) -I$(top_builddir)/psi \
		-I$(top_srcdir) -I$(top_srcdir)/psi

if WITH_READLINE
LIBADD_READLINE = -lreadline
else
LIBADD_READLINE =
endif

PSILIBRARIES	= $(top_builddir)/psi/libpsi.a \
		$(top_builddir)/replace/libreplace.a $(LIBADD_READLINE) -lgsl -lgslcblas -lm 

noinst_HEADERS	= psi.h

bin_PROGRAMS	= $(TESTS)

psi_SOURCES	= psi.c
psi_LDADD	= $(PSILIBRARIES)

seq_SOURCES	= seq.c
seq_LDADD	= $(PSILIBRARIES)

grid_SOURCES	= grid.c
grid_LDADD	= $(PSILIBRARIES)

energy_SOURCES	= energy.c
energy_LDADD	= $(PSILIBRARIES)

rng_SOURCES	= rng.c
rng_LDADD	= $(PSILIBRARIES)

gibbs_SOURCES	= gibbs.c
gibbs_LDADD	= $(PSILIBRARIES)

mcmc_SOURCES	= mcmc.c
mcmc_LDADD	= $(PSILIBRARIES)

ds_SOURCES	= ds.c
ds_LDADD	= $(PSILIBRARIES)

bf_SOURCES	= bf.c
bf_LDADD	= $(PSILIBRARIES)

mc_SOURCES	= mc.c
mc_LDADD	= $(PSILIBRARIES)

nsv_SOURCES	= nsv.c
nsv_LDADD	= $(PSILIBRARIES)

sim_SOURCES	= sim.c
sim_LDADD	= $(PSILIBRARIES)

