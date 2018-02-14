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
gmpmee_spowm_init(gmpmee_spowm_tab table, size_t len, mpz_t modulus,
		  size_t block_width)
{
  size_t i, j;
  size_t tab_len;  /* Size of a subtable. */
  mpz_t *t;        /* Temporary variable for subtable. */

  table->len = len;
  table->block_width = block_width;
  if (len < block_width) {
    table->block_width = len;
  }
  table->tabs_len = (len + block_width - 1) / block_width;

  mpz_init(table->modulus);
  mpz_set(table->modulus, modulus);

  /* Allocate and initialize space for pointers to tables. */
  table->tabs = (mpz_t **)malloc(table->tabs_len * sizeof(mpz_t *));

  tab_len = 1 << block_width;
  for (i = 0; i < table->tabs_len; i++)
    {

      /* Last block may be more narrow than the other, but they are
         never zero. */
      if (i == table->tabs_len - 1
	  && len - (table->tabs_len - 1) * block_width < block_width)
	{
	  block_width = len - (table->tabs_len - 1) * block_width;
	  tab_len = 1 << block_width;
	}

      /* Allocate and initialize a table. */
      table->tabs[i] = (mpz_t *)malloc(tab_len * sizeof(mpz_t));
      t = table->tabs[i];

      /* Initialize mpz_t's. */
      for (j = 0; j < tab_len; j++)
	{
	  mpz_init(t[j]);
	}
    }
}
