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
gmpmee_millerrabin_next_cand(gmpmee_millerrabin_state state)
{

  /* Add two if n is odd. */
  if (mpz_tstbit(state->n, 0))
    {
      mpz_add_ui(state->n, state->n, 2L);
    }
  else /* Add one otherwise. */
    {
      mpz_add_ui(state->n, state->n, 1L);
    }

  /* Add two until trial divisions do not reveal the number as
     composite. */
  while (gmpmee_millerrabin_trial(state->n) == 0)
    {
      mpz_add_ui(state->n, state->n, 2L);
    }

  /* Update the state and define q and k such that n = q*2^k+1. */
  mpz_sub_ui(state->n_minus_1, state->n, 1L);
  state->k = mpz_scan1(state->n_minus_1, 0L);
  mpz_tdiv_q_2exp(state->q, state->n_minus_1, state->k);
}
