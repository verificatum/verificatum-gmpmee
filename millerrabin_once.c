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
gmpmee_millerrabin_once(gmpmee_millerrabin_state state, mpz_t base)
{
  unsigned long int i;

  if (mpz_cmp_ui(state->n, 4) < 0) {
    if (mpz_cmp_ui(state->n, 1) > 0) {
      return 1;
    } else {
      return 0;
    }
  }

  mpz_powm(state->y, base, state->q, state->n);

  if (mpz_cmp_ui(state->y, 1L) == 0 || mpz_cmp(state->y, state->n_minus_1) == 0)
    {
      return 1;
    }

  for (i = 1; i < state->k; i++)
    {

      mpz_powm_ui(state->y, state->y, 2L, state->n);

      if (mpz_cmp_ui(state->y, 1L) == 0)
	{
	  return 0;
	}

      if (mpz_cmp(state->y, state->n_minus_1) == 0)
	{
	  return 1;
	}
    }
  return 0;
}
