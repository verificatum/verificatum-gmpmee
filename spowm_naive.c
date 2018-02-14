/*
 * Copyright 2008 2009 2010 2011 2013 2014 2015 2016 Douglas Wikstrom
 *
 * This file is part of GMP Modular Exponentiation Extension (GMPMEE).
 *
 * GMPMEE is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GMPMEE is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GMPMEE. If not, see <http://www.gnu.org/licenses/>.
 */

#include <gmp.h>
#include "gmpmee.h"

void
gmpmee_spowm_naive(mpz_t rop, mpz_t *bases, mpz_t *exponents, size_t len,
		   mpz_t modulus)
{
  size_t i;
  mpz_t tmp;

  mpz_set_ui(rop, 1);
  mpz_init(tmp);

  for (i = 0; i < len; i++)
    {
      mpz_powm(tmp, bases[i], exponents[i], modulus);
      mpz_mul(rop, rop, tmp);
      mpz_mod(rop, rop, modulus);
    }

  mpz_clear(tmp);
}
