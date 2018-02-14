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

void
mpz_probab_safe_prime_p_next(mpz_t rop, mpz_t n, int reps)
{
  int increased = 0;

  mpz_set(rop, n);

  if (mpz_cmp_ui(rop, 5) < 0)
    {
      mpz_set_ui(rop, 5);
    }
  else if (mpz_cmp_ui(rop, 7) < 0)
    {
      mpz_set_ui(rop, 7);
    }
  else
    {

      /* Make sure that rop is odd. */
      if (!mpz_tstbit(rop, 0))
        {
          mpz_add_ui(rop, rop, 1L);
          increased = 1;
        }

      /* Make sure that m is odd, where rop=2m+1. */
      if (!mpz_tstbit(rop, 1))
        {
          mpz_add_ui(rop, rop, 2L);
          increased = 1;
        }

      /* If both rop and m were already odd, then we add 4. */
      if (!increased)
        {
          mpz_add_ui(rop, rop, 4L);
        }

      while (mpz_probab_safe_prime_p(rop, reps) == 0)
        {
          mpz_add_ui(rop, rop, 4L);
        }
    }
}
