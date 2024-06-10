# GMP Modular Exponentiation Extension (GMPMEE)

## Overview

This is a minor extension of the Gnu Multiprecision Library (GNU MP).
It adds simultaneous modular exponentiation and fixed base modular
exponentiation functionality to the set of integer functions (the
mpz-functions), as well as special purpose primality testing routines.

GMP does contain primality testing routines, but these do not use
cryptographically strong randomness and they do not allow fast testing
or sieving for safe-primality.

Note that **no attempt is made to make this secure against
side-channel attacks**, since there is no need for this type of
protection in the applications for which this library was implemented.

For a detailed account of such algorithms, a good source is [Handbook
of Applied Cryptography](http://www.cacr.math.uwaterloo.ca/hac),
*Menezes, Oorshot, and Vanstone*, which is available for free.

The following assumes that you are using a distribution. Developers
should also read `README_DEV.md`.


## Building

If you are using a distribution, then you simply use

        ./configure
        make

to build the library and an executable `gmpmee` that allows testing
and benchmarking some of the routines.

If you prefer to use the [Clang compiler](https://clang.llvm.org) in
place of GCC for the native code, then you may use `./configure
CC=clang` instead of the above to enable it.

**Caution: Please understand that although it seems that Clang works
as well as GCC, switching compiler is a large change for mature
software.**


## Installing

Use

        make install

to install the library `libgmpmee.{la,a,so}` and the corresponding
header file `gmpmee.h` in the standard locations. See `INSTALL` for
details on other ways to invoke `./configure`, e.g., to use a
user-local installation.

You may need to run `sudo /sbin/ldconfig` on some platforms which have
flawed implementations of the cache that stores locations of
libraries.

## Usage

If you have done a standard install you may use the library by
including `gmpmee.h` and adding the flags `-lgmp -lgmpmee`. We could
for example compile a program `foo.c` using the library as

        gcc foo.c -o foo -lgmpmee -lgmp

If you have done a non-standard installation you may need to update
some environment variables.


## API Documentation

You may use

        make api

to build also some documentation using Doxygen (this assumes you have
installed `doxygen`). The API is not installed anywhere. You can copy
it to any location.


## Reporting Bugs

Minor bugs should be reported in the repository system as issues or
bugs. Security critical bugs, vulnerabilities, etc should be reported
directly to the
[https://www.verificatum.org](https://www.verificatum.org). We will
make best effort to disclose the information in a responsible way
before the finder gets proper credit.
