/*
 * Copyright 2008 2009 2010 2011 2012 Torbjorn Granlund, Douglas Wikstrom
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

int
mpz_probab_safe_prime_p(mpz_t n, int reps)
{
  mpz_t m;
  int res;

  if (!mpz_probab_prime_p(n, reps + 1))
    {
      return 0;
    }

  mpz_init(m);
  mpz_sub_ui(m, n, 1L);
  mpz_div_ui(m, m, 2L);

  res = mpz_probab_prime_p(m, reps + 1);

  mpz_clear(m);

  return res;
}
