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

#include <gmp.h>
#include "gmpmee.h"

void
gmpmee_fpowm_precomp(gmpmee_fpowm_tab table, mpz_t basis)
{
  size_t i;
  size_t block_width = table->spowm_table->block_width;
  mpz_t *bases = gmpmee_array_alloc_init(block_width);

  mpz_t eb;

  /* eb = 2^(table->stretch) */
  mpz_init(eb);
  mpz_setbit(eb, table->stretch);

  mpz_set(bases[0], basis);
  for (i = 1; i < block_width; i++) {
    mpz_powm(bases[i], bases[i - 1], eb, table->spowm_table->modulus);
  }

  gmpmee_spowm_precomp(table->spowm_table, bases);

  mpz_clear(eb);
  gmpmee_array_clear_dealloc(bases, block_width);
}
