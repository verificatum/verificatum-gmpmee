
# Copyright 2008 2009 2010 2011 2013 2014 2015 2016 Douglas Wikstrom
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

AC_PREREQ([2.69])
AC_INIT([verificatum-gmpmee],[m4_esyscmd([m4 .version.m4])],[info@verificatum.com])
AC_LANG_PUSH([C])
AM_INIT_AUTOMAKE([foreign tar-ustar -Wall -Werror -Wno-portability])
AC_CONFIG_SRCDIR([spowm_init.c])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# Disable any default flags set by AC_PROG_CC. This still allows the
# user to specify CFLAGS.
CFLAGSTEST=`echo ${CFLAGS} | tr -d " "`
if test x${CFLAGSTEST} = x;
then
   # This is needed since AC_PROG_CC is implemented incorrectly. It
   # checks if CFLAGS is undefined or defined to the empty string. It
   # should check if CFLAGS is white space.
   CFLAGS="-fmessage-length=0"
fi 

# Checks for programs.
AC_PROG_CC

# We are building a library.
LT_INIT

# Checks for libraries.
AC_CHECK_LIB(gmp, __gmpz_init, ,
       [AC_MSG_ERROR(["GNU MP library not found, see http://gmplib.org"])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
AC_CHECK_HEADERS([gmp.h], ,
       [AC_MSG_ERROR(["GNU MP header not found, see http://gmplib.org"])])

# Compile a small program that extracts the compiler flags used to
# compile GMP.
${CC} extract_GMP_CFLAGS.c -o extract_GMP_CFLAGS

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
