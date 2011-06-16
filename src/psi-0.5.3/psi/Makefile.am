## Makefile.am -- Process this file with automake to produce Makefile.in
## Copyright (C) 2006 Sang Chul Choi
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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
## MA 02110-1301, USA.

MAINTAINERCLEANFILES	= Makefile.in

INCLUDES		= -I$(top_builddir) -I$(top_srcdir)

pkginclude_HEADERS 	= common.h error.h psi.h seq.h grid.h energy.h \
                          rng.h gslwrap.h gibbs.h mcmc.h io.h bf.h ds.h \
                          mc.h nsv.h sim.h defs.h

lib_LIBRARIES 		= libpsi.a
libpsi_a_SOURCES 	= error.c psi.c seq.c grid.c energy.c rng.c \
			  gslwrap.c gibbs.c mcmc.c io.c bf.c ds.c \
			  mc.c nsv.c sim.c \
			  xmalloc.c xstrdup.c xstrerror.c

