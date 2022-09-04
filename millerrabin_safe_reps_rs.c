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
gmpmee_millerrabin_safe_reps_rs(gmp_randstate_t rstate,
				gmpmee_millerrabin_safe_state state,
				int reps)
{
  int i;
  int res = 1;
  mpz_t nn_minus_1;
  mpz_t nbase;
  mpz_t mn_minus_1;
  mpz_t mbase;

  mpz_init(nn_minus_1);
  mpz_sub_ui(nn_minus_1, state->nstate->n, 1);
  mpz_init(nbase);

  mpz_init(mn_minus_1);
  mpz_sub_ui(mn_minus_1, state->mstate->n, 1);
  mpz_init(mbase);

  /* Repeat an additional time. The union bound then gives the
     expected provable error probability. */
  for (i = 0; i < reps + 1; i++) {

    /* Random base in [2,n-2] */
    mpz_urandomm(nbase, rstate, nn_minus_1);
    if (mpz_cmp_ui(nbase, 2) < 0)
      {
        mpz_set_ui(nbase, 2);
      }

    /* FIXME: GCC + libtool is currently broken. A fixed-size array
       within a struct is valid C code and sizeof() should determine
       the correct size of the resulting struct, including any
       padding, but without disabling this libtool gives a warning. It
       would still be better to rewrite the code to avoid the buggy
       error. This is however quite difficult to do in a backwards
       compatible way. */

#ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
    if (gmpmee_millerrabin_once(state->nstate, nbase) == 0)
      {
	res = 0;
	break;
      }
#ifndef __clang__
#pragma GCC diagnostic pop
#endif

    /* Random base in [2,m-2] */
    mpz_urandomm(mbase, rstate, mn_minus_1);
    if (mpz_cmp_ui(mbase, 2) < 0)
      {
        mpz_set_ui(mbase, 2);
      }

#ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
    if (gmpmee_millerrabin_once(state->mstate, mbase) == 0)
      {
	res = 0;
	break;
      }
#ifndef __clang__
#pragma GCC diagnostic pop
#endif
  }

  mpz_clear(nn_minus_1);
  mpz_clear(nbase);

  mpz_clear(mn_minus_1);
  mpz_clear(mbase);

  return res;
}
