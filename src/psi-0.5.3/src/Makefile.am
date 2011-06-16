## Makefile.am -- Process this file with automake to produce Makefile.in
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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
## MA 02110-1301, USA.

MAINTAINERCLEANFILES = Makefile.in

INCLUDES	= -I$(top_builddir) -I$(top_builddir)/psi \
		-I$(top_srcdir) -I$(top_srcdir)/psi

if WITH_READLINE
LIBADD_READLINE = -lreadline
else
LIBADD_READLINE =
endif

noinst_HEADERS	= psi.h defs.h

bin_PROGRAMS	= psi
psi_SOURCES	= psi.c ## psi_builtin.c psi_repl.c psi_syntax.c
psi_LDADD	= $(top_builddir)/psi/libpsi.a \
		$(top_builddir)/replace/libreplace.a $(LIBADD_READLINE) -lgsl -lgslcblas -lm 

