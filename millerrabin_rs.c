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
gmpmee_millerrabin_rs(gmp_randstate_t rstate, mpz_t n, int reps)
{

  int res;
  gmpmee_millerrabin_state state;

  if (mpz_cmp_ui(n, 4) < 0)
    {
      return mpz_tstbit(n, 1);
    }
  else if (gmpmee_millerrabin_trial(n) == 0)
    {
      return 0;
    }
  else
    {
      gmpmee_millerrabin_init(state, n);
      res = gmpmee_millerrabin_reps_rs(rstate, state, reps);
      gmpmee_millerrabin_clear(state);
      return res;
    }
}
