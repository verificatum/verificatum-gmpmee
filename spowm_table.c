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

/*
 * Returns the index'th bit of each of the first block_width integers
 * in the array. The least significant bit in the output is the bit
 * extracted from the first integer in the input array.
 */
static int
getbits(mpz_t *op, int index, size_t block_width)
{
  int i;
  int bits = 0;

  for (i = block_width - 1; i >= 0; i--)
    {
      bits <<= 1;
      if (mpz_tstbit(op[i], index))
	{
	  bits |= 1;
	}
    }
  return bits;
}

void
gmpmee_spowm_table(mpz_t rop, gmpmee_spowm_tab table, mpz_t *exponents)
{
  size_t i;
  int index;
  int mask;
  size_t bitlen;
  mpz_t *exps;
  size_t max_exponent_bitlen;
  size_t len = table->len;
  size_t tabs_len = table->tabs_len;
  size_t block_width = table->block_width;
  size_t last_block_width = len - (tabs_len - 1) * block_width;
  mpz_t **tabs = table->tabs;

  /* Compute the maximal bit length among the exponents. */
  max_exponent_bitlen = 0;
  for (i = 0; i < len; i++)
    {
      bitlen = mpz_sizeinbase(exponents[i], 2);
      if (bitlen > max_exponent_bitlen)
	{
	  max_exponent_bitlen = bitlen;
	}
    }

  /* Initialize result variable. */
  mpz_set_ui(rop, 1);

  /* Execute simultaneous square-and-multiply. */
  for (index = max_exponent_bitlen - 1; index >= 0; index--)
    {

      /* Square ... */
      mpz_mul(rop, rop, rop);
      mpz_mod(rop, rop, table->modulus);

      /* ... and multiply. */
      i = 0;
      exps = exponents;
      while (i < tabs_len)
	{
	  if (i == tabs_len - 1)
	    {
	      mask = getbits(exps, index, last_block_width);
	    }
	  else
	    {
              /* This is never executed if there is a single table. */
              /* LCOV_EXCL_START */
	      mask = getbits(exps, index, block_width);
              /* LCOV_EXCL_STOP */
	    }

	  mpz_mul(rop, rop, tabs[i][mask]);
          mpz_mod(rop, rop, table->modulus);
	  i++;
	  exps += block_width;
	}
    }
}
