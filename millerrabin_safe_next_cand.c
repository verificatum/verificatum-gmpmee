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
gmpmee_millerrabin_safe_next_cand(gmpmee_millerrabin_safe_state state)
{
  int increased = 0;

  /* Make sure that n is odd. */
  if (!mpz_tstbit(state->nstate->n, 0))
    {
      mpz_add_ui(state->nstate->n, state->nstate->n, 1L);
      increased = 1;
    }

  /* Make sure that m is odd, where n=2m+1. */
  if (!mpz_tstbit(state->nstate->n, 1))
    {
      mpz_add_ui(state->nstate->n, state->nstate->n, 2L);
      increased = 1;
    }

  /* If both n and m were already odd, then we add 4. */
  if (!increased) {
    mpz_add_ui(state->nstate->n, state->nstate->n, 4L);
  }

  /* Add four until trial divisions do not reveal the number as
     non safe-prime. */
  while (gmpmee_millerrabin_safe_trial(state->nstate->n) == 0)
    {
      mpz_add_ui(state->nstate->n, state->nstate->n, 4L);
    }

  /* Update the state for testing of n, and define q and k such that
     n=q*2^k+1. (q and k are local to the state) */
  mpz_sub_ui(state->nstate->n_minus_1, state->nstate->n, 1L);
  state->nstate->k = mpz_scan1(state->nstate->n_minus_1, 0L);
  mpz_tdiv_q_2exp(state->nstate->q, state->nstate->n_minus_1, state->nstate->k);

  /* Update the state for testing of m, where n=2m+1, and define q and
     k such that m=q*2^k+1. (q and k are local to the state) */
  mpz_div_ui(state->mstate->n, state->nstate->n_minus_1, 2);
  mpz_sub_ui(state->mstate->n_minus_1, state->mstate->n, 1L);
  state->mstate->k = mpz_scan1(state->mstate->n_minus_1, 0L);
  mpz_tdiv_q_2exp(state->mstate->q, state->mstate->n_minus_1, state->mstate->k);
}
