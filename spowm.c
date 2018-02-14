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
 * Computes a "theoretical" optimal block width for a given exponent
 * length. This is typically not optimal in practice, so we only use
 * this for values outside the table below.
 */
static int
theoretical_block_width(int exponents_bitlen) {
  int res = -1;
  int est;
  int w;
  int opt_est;
  int opt_w = 1;

  opt_est = exponents_bitlen;

  for (w = 2; w < 50; w++)
    {

      /* Precomputation: 2^w - w - 1
       * Computation: exponents_bitlen/w */
      est = (1 << w) - w - 1 + exponents_bitlen/w;

      if (est > opt_est)
	{
	  res = opt_w;
          break;
	}
      else
	{
	  opt_est = est;
	  opt_w = w;
	}
    }
  return res;
}

#define GMPMEE_SPOWM_ROWS 7
#define GMPMEE_SPOWM_COLUMNS 8

/* We have separate tables for the following bitlengths. */
size_t modulus_bitlens[] = {64, 128, 256, 512, 1024, 2048, 4096};

/* Table of "practical optima" for some block widths. Zero is used
   when no value is set. This should really be tuned correctly for
   each individual platform. */
size_t best_block_widths[GMPMEE_SPOWM_ROWS][GMPMEE_SPOWM_COLUMNS] =
{
/*            5    6    7     8     9    10    11  12  */
/*   64 */ {100, 150, 350, 1100, 1100,    0,    0,  0},
/*  128 */ {100, 150, 350, 1000, 1350,    0,    0,  0},
/*  256 */ {100, 150, 450, 4450,    0,    0,    0,  0},
/*  512 */ {100, 200, 500, 1700, 5000,    0,    0,  0},
/* 1024 */ {100, 100, 500, 1000, 2500, 6000,    0,  0},
/* 2048 */ {100, 150, 450, 1000, 2000, 4500, 8200,  0},
/* 4096 */ {100, 200, 350,  900, 2000, 4400, 7300,  0}
};

void
gmpmee_spowm(mpz_t rop, mpz_t *bases, mpz_t *exponents, size_t len,
	     mpz_t modulus)
{
  size_t i;
  int row;
  int column;
  size_t block_width;
  size_t modulus_bitlen;
  size_t bitlen;
  size_t max_exponent_bitlen;

  /* Process in a single batch. */
  size_t batch_len = len;

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

  /* Find bit length of modulus. */
  modulus_bitlen = mpz_sizeinbase(modulus, 2);

  /* Lookup block-width table for the given modulus bit length. */
  for (row = 0; row < GMPMEE_SPOWM_ROWS; row++)
    {
      if (modulus_bitlens[row] > modulus_bitlen)
	{
	  break;
	}
    }
  if (row > 0)
    {
      row--;
    }

  /* Lookup "optimal" block-width. */
  for (column = 0; column < GMPMEE_SPOWM_COLUMNS; column++)
    {
      if (best_block_widths[row][column] > max_exponent_bitlen)
	{
	  break;
	}
    }
  if (column == GMPMEE_SPOWM_COLUMNS)
    {
      block_width = theoretical_block_width(max_exponent_bitlen) + 2;
    }
  else
    {
      if (column > 0)
	{
	  column--;
	}
      block_width = column + 5;
    }

  gmpmee_spowm_block_batch(rop, bases, exponents, len, modulus, block_width,
			   batch_len);
}
