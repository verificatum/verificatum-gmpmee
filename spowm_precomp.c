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

#include <stdlib.h>
#include <gmp.h>
#include "gmpmee.h"

void
gmpmee_spowm_precomp(gmpmee_spowm_tab table, mpz_t *bases)
{
  size_t i, j;
  size_t tabs_len = table->tabs_len;
  size_t block_width = table->block_width;
  int mask;
  int one_mask;
  mpz_t *t;

  for (i = 0; i < tabs_len; i++)
    {
      /* Last block may have smaller width, but it is never zero. */
      if (i == tabs_len - 1)
        {
          block_width = table->len - (tabs_len - 1) * block_width;
        }

      /* Current subtable. */
      t = table->tabs[i];

      /* Initialize current subtable with all trivial products. */
      mpz_set_ui(t[0], 1);

      mask = 1;
      for (j = 0; j < block_width; j++)
        {
          mpz_set(t[mask], bases[j]);
          mask <<= 1;
        }

      /* Initialize current subtable with all non-trivial products. */
      for (mask = 1; mask < (1 << block_width); mask++)
        {
          one_mask = mask & (-mask);
          mpz_mul(t[mask], t[mask ^ one_mask], t[one_mask]);
          mpz_mod(t[mask], t[mask], table->modulus);
        }

      bases += block_width;
    }
}
