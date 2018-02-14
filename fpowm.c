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
 * Let op = (x_0,..,x_t), where x_i is a binary string with stretch
 * number of bits (except that x_t may have more bits) and
 * t=block_width. This function returns (x_{0,i},...,x_{t,i}), where
 * x_{j,i} is the ith bit of x_j. Each x_j except x_t is considered as
 * padded with zeros if i > stretch is asked for.
 *
 * @param op Integer from which bits are derived.
 * @param index Bit index in each sequence of stretch (or more) bits.
 * @param block_width Number of integers considered in parallel.
 * @param stretch Number of bits in each integer.
 */
static int
getbits(mpz_t op, int index, size_t block_width, size_t stretch)
{
  int i;
  int bits = 0;

  if (((size_t) index) < stretch)
    {
      for (i = block_width - 1; i >= 0; i--)
	{
	  bits <<= 1;
	  if (mpz_tstbit(op, i * stretch + index))
	    {
	      bits |= 1;
	    }
	}
      return bits;

    }
  else
    {
      if (mpz_tstbit(op, (block_width - 1) * stretch + index)) {
	bits |= 1;
	bits <<= (block_width - 1);
      }
      return bits;

    }
}

void
gmpmee_fpowm(mpz_t rop, gmpmee_fpowm_tab table, mpz_t exponent)
{
  int index;
  int mask;
  size_t max_exponent_bitlen;
  size_t block_width = table->spowm_table->block_width;
  mpz_t **tabs = table->spowm_table->tabs;


  /* Compute the maximal bit length among the exponents. */
  max_exponent_bitlen = mpz_sizeinbase(exponent, 2);
  if (max_exponent_bitlen < block_width * table->stretch) {
    max_exponent_bitlen = table->stretch;
  } else {
    max_exponent_bitlen -= (block_width - 1) * table->stretch;
  }


  mpz_set_ui(rop, 1);

  /* Execute square-and-multiply. */
  for (index = max_exponent_bitlen - 1; index >= 0; index--)
    {

      /* Square ... */
      mpz_mul(rop, rop, rop);
      mpz_mod(rop, rop, table->spowm_table->modulus);

      /* and multiply */
      mask = getbits(exponent, index, block_width, table->stretch);
      mpz_mul(rop, rop, tabs[0][mask]);
      mpz_mod(rop, rop, table->spowm_table->modulus);

    }
}
