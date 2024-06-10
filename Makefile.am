
# Copyright 2008 2009 2010 2011 2013 Torbjorn Granlund, Douglas Wikstrom
#
# This file is part of GMP Modular Exponentiation Extension (GMPMEE).
#
# GMPMEE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GMPMEE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GMPMEE. If not, see <http://www.gnu.org/licenses/>.

.PHONY: clean cleanapi

BINDIR = bin

ACLOCAL_AMFLAGS=-I m4

# Figuring out the best compiler/building flags to use for a given
# platform is difficult. This project uses GMP to implement large
# integers, so the flags used by GMP to compile should be a good
# default on all our target platforms.
#
# The authors of GMP had the good taste to provide these flags as the
# macro __GMP_CFLAGS, so writing a program to print the flags is
# trivial. Our program is compiled during configuration.

GMP_CFLAGS := $(shell ./extract_GMP_CFLAGS)

OPTIONAL_FLAGS=-O3

# Use this instead to run coverage analysis.
# OPTIONAL_FLAGS=-fprofile-arcs -ftest-coverage

# We use pedantic flags, strip the optimization flag of GMP, and
# insert our own level of optimization.
AM_CFLAGS := -Wall -W -Werror $(shell echo ${GMP_CFLAGS} | sed -e "s/-O[O12345]//") $(OPTIONAL_FLAGS)


AM_LDFLAGS = -lgmp

scriptmacros.m4:
	@printf "define(M4_VERSION, $(VERSION))dnl\n" > scriptmacros.m4

$(BINDIR)/gmpmee-info: scriptmacros.m4 gmpmee-info.src
	mkdir -p $(BINDIR)
	cat scriptmacros.m4 gmpmee-info.src | m4 > $(BINDIR)/gmpmee-info
	chmod +x $(BINDIR)/gmpmee-info


# Actual target.
lib_LTLIBRARIES = libgmpmee.la
libgmpmee_la_SOURCES = spowm.c spowm_clear.c spowm_precomp.c spowm_init.c spowm_table.c spowm_block_batch.c spowm_naive.c array_alloc.c array_clear_dealloc.c array_urandomb.c array_alloc_init.c fpowm.c fpowm_clear.c fpowm_init.c fpowm_precomp.c fpowm_init_precomp.c millerrabin_init.c millerrabin_next_cand.c millerrabin_clear.c millerrabin_trial.c millerrabin_once.c millerrabin_reps_rs.c millerrabin_rs.c millerrabin_next_rs.c millerrabin_safe_init.c millerrabin_safe_next_cand.c millerrabin_safe_clear.c millerrabin_safe_trial.c millerrabin_safe_reps_rs.c millerrabin_safe_rs.c millerrabin_safe_next_rs.c probab_prime_p_next.c probab_safe_prime_p.c probab_safe_prime_p_next.c

libgmpmee_la_LIBADD = -lgmp
gmpmee_LDADD = libgmpmee.la

include_HEADERS = gmpmee.h
gmpmee_SOURCES = gmpmee.c gmpmee.h
bin_PROGRAMS = gmpmee
dist_bin = $(BINDIR)/gmpmee-info
dist_bin_SCRIPTS = $(BINDIR)/gmpmee-info

dist_noinst_DATA = extract_GMP_CFLAGS.c doxygen.cfg trialdiv_32.c trialdiv_64.c trialdiv_safe_32.c trialdiv_safe_64.c .version.m4 gmpmee-info.src

all-local: check_info.stamp

check_info.stamp: gmpmee
	@echo ""
	@echo "================================================================"
	@echo ""
	@echo "          INSTALLING AND VERIFYING INSTALLATION"
	@echo ""
	@echo " To install, you now use the standard 'make install' (or rather"
	@echo " 'sudo make install' on most Linuxes), but before you do that "
	@echo " it is prudent to use 'make check' to verify your installation."
	@echo ""
	@echo "================================================================"
	@echo ""
	@touch check_info.stamp

install-data-hook:
	@echo ""
	@echo "================================================================"
	@echo ""
	@echo "                      TESTING"
	@echo ""
	@echo " Please consult the README file for information about how to"
	@echo " use the gmpmee executable to test the library, or simply "
	@echo " execute 'gmpmee' to print usage information."
	@echo ""
	@echo "================================================================"
	@echo ""

check:
	@./gmpmee 800

# Uncomment OPTIONAL_FLAGS above to enable instrumentation.
coverage: all
	@./gmpmee 1000
	@rm -rf coverage
	@mkdir -p coverage
	@lcov --capture --directory . --output-file coverage/coverage.info
	@genhtml coverage/coverage.info --output-directory coverage/html

api:
	@mkdir -p html_api
	doxygen doxygen.cfg
	mv html_api/html api
	rm -rf html_api

cleancoverage:
	@rm -rf *.gcno *.gcda coverage coverage.info

cleanapi:
	rm -rf api

clean-local: cleanapi cleancoverage
	find . -name "*~" -delete
	rm -rf *.stamp scriptmacros.m4 $(BINDIR)
