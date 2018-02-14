/*
 * Copyright 2008 2009 2010 2011 2013 Torbjorn Granlund, Douglas Wikstrom
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
mpz_probab_prime_p_next(mpz_t rop, mpz_t n, int reps)
{
  mpz_set(rop, n);

  if (mpz_cmp_ui(rop, 2) < 0)
    {
      mpz_set_ui(rop, 2);
    }
  else
    {
      /* Add two if n is odd. */
      if (mpz_tstbit(rop, 0))
        {
          mpz_add_ui(rop, rop, 2L);
        }
      else /* Add one otherwise. */
        {
          mpz_add_ui(rop, rop, 1L);
        }

      while (mpz_probab_prime_p(rop, reps) == 0)
        {
          mpz_add_ui(rop, rop, 2L);
        }
    }
}
