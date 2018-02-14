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

int
gmpmee_millerrabin_reps_rs(gmp_randstate_t rstate,
			   gmpmee_millerrabin_state state,
			   int reps)
{
  int i;
  int res;
  mpz_t base;
  mpz_t n_minus_1;

  mpz_init(base);
  mpz_init(n_minus_1);
  mpz_sub_ui(n_minus_1, state->n, 1);

  res = 1;
  for (i = 0; res && i < reps; i++) {

    /* Almost random base in [2,n-2] */
    mpz_urandomm(base, rstate, n_minus_1);
    if (mpz_cmp_ui(base, 2) < 0)
      {
        mpz_set_ui(base, 2);
      }

    res = gmpmee_millerrabin_once(state, base);
  }

  mpz_clear(base);
  mpz_clear(n_minus_1);

  return res;
}
