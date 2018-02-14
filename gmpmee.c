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

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gmp.h>
#include "gmpmee.h"

int
gmpmee_done(long start_time, long interval)
{
  return (clock() - start_time > interval * (CLOCKS_PER_SEC / 1000));
}

void
test_miller_rabin_n(int call, gmp_randstate_t rstate, mpz_t n)
{
  int res;
  int reps = 20;
  mpz_t rop;
  mpz_t gmprop;

  mpz_init(rop);
  mpz_init(gmprop);

  switch (call)
    {
    case 0:
      res = mpz_probab_prime_p(n, reps) ? 1 : 0;
      assert (res == gmpmee_millerrabin_rs(rstate, n, reps));
      break;

    case 1:
      res = mpz_probab_safe_prime_p(n, reps) ? 1 : 0;
      assert (res == gmpmee_millerrabin_safe_rs(rstate, n, reps));
      break;

    case 2:
      mpz_probab_prime_p_next(gmprop, n, reps);
      gmpmee_millerrabin_next_rs(rop, rstate, n, reps);
      assert (mpz_cmp(gmprop, rop) == 0);
      break;

    default:
      mpz_probab_safe_prime_p_next(gmprop, n, reps);
      gmpmee_millerrabin_safe_next_rs(rop, rstate, n, reps);
      assert (mpz_cmp(gmprop, rop) == 0);
    }

  mpz_clear(gmprop);
  mpz_clear(rop);
}

void
test_miller_rabin(int call, long test_time)
{
  int t;
  int i;
  gmp_randstate_t rstate;
  mpz_t n;
  int bit_length = 256;

  gmp_randinit_default(rstate);
  mpz_init(n);
  t = clock();

  /* Check all small numbers for sanity. */
  for (i = 0; i < 20000; i++) {
    mpz_set_ui(n, i);
    test_miller_rabin_n(call, rstate, n);
  }

  mpz_urandomb(n, rstate, bit_length);
  do
    {
      mpz_add_ui(n, n, 1L);
      test_miller_rabin_n(call, rstate, n);
    }
  while (!gmpmee_done(t, test_time));

  mpz_clear(n);
  gmp_randclear(rstate);
}

void
test_spowm_modulus_bitlen(int modulus_bitlen)
{

  int i;
  int len;
  int exponents_bitlen;

  gmp_randstate_t state;
  mpz_t modulus;
  mpz_t *bases;
  mpz_t *exponents;

  mpz_t naive_res;
  mpz_t spowm_res;

  /* Initialize random state. */
  gmp_randinit_default(state);

  mpz_init(modulus);
  mpz_init(naive_res);
  mpz_init(spowm_res);

  len = 1;

  /* Generate modulus. */
  do
    {
      mpz_urandomb(modulus, state, modulus_bitlen);
    }
  while (mpz_cmp_ui(modulus, 0) == 0);


  /* Allocate space. */
  bases = gmpmee_array_alloc_init(len);
  exponents = gmpmee_array_alloc_init(len);


  /* Generate bases. */
  gmpmee_array_urandomb(bases, len, state, modulus_bitlen);
  for (i = 0; i < len; i++)
    {
      mpz_mod(bases[i], bases[i], modulus);
    }

  /* Generate exponents for a few bit lengths and perform the
     test. */
  exponents_bitlen = modulus_bitlen / 2;
  if (exponents_bitlen == 0)
    {
      exponents_bitlen = 1;
    }

  do
    {

      /* Generate exponents. */
      gmpmee_array_urandomb(exponents, len, state, exponents_bitlen);

      /* Compute in both ways. */
      gmpmee_spowm_naive(naive_res, bases, exponents, len, modulus);
      gmpmee_spowm(spowm_res, bases, exponents, len, modulus);

      /* Compare results. */
      assert(mpz_cmp(spowm_res, naive_res) == 0);

      /* Increase bit lengths of exponents. */
      exponents_bitlen <<= 1;
    }
  while (exponents_bitlen < 3 * modulus_bitlen);

  gmpmee_array_clear_dealloc(exponents, len);
  gmpmee_array_clear_dealloc(bases, len);

  /* Make sure corner case is handled. */
  len = 2 * len + 1;
}

void
test_spowm(long test_time)
{
  int t;
  int i;

  t = clock();

  do
    {
      for (i = 1; i < 10; i++)
        {
          test_spowm_modulus_bitlen(i);
        }
      test_spowm_modulus_bitlen(1024);
      test_spowm_modulus_bitlen(10000);

    } while (!gmpmee_done(t, test_time));
}

void
test_fpowm(long test_time)
{

  int t;
  int modulus_bitlen = 1024;
  int exponent_bitlen;
  int block_width = 0;

  gmp_randstate_t state;
  mpz_t modulus;
  mpz_t basis;
  mpz_t exponent;

  mpz_t naive_res;
  mpz_t fpowm_res;
  gmpmee_fpowm_tab table;

  /* Initialize random state. */
  gmp_randinit_default(state);

  mpz_init(modulus);
  mpz_init(basis);
  mpz_init(exponent);
  mpz_init(naive_res);
  mpz_init(fpowm_res);

  exponent_bitlen = 1;

  t = clock();

  do
    {

      if (block_width == 0) {
        block_width++;
      }

      /* Generate modulus. */
      do
	{
	  mpz_urandomb(modulus, state, modulus_bitlen);
	}
      while (mpz_cmp_ui(modulus, 0) == 0);

      /* Generate basis. */
      do
	{
	  mpz_urandomb(basis, state, modulus_bitlen);
	  mpz_mod(basis, basis, modulus);
	}
      while (mpz_cmp_ui(basis, 1) <= 0);

      gmpmee_fpowm_init_precomp(table,
                                basis,
                                modulus,
                                block_width,
                                exponent_bitlen);
      do
	{

	  /* Generate exponent. */
	  mpz_urandomb(exponent, state, exponent_bitlen);

	  /* Compute in both ways. */
	  gmpmee_fpowm(fpowm_res, table, exponent);
	  mpz_powm(naive_res, basis, exponent, modulus);

	  /* Compare results. */
	  assert(mpz_cmp(fpowm_res, naive_res) == 0);

	  /* Increase bit length of exponent. */
	  exponent_bitlen <<= 1;
	}
      while (exponent_bitlen < 3 * modulus_bitlen);

      gmpmee_fpowm_clear(table);

      block_width = (block_width + 1) % modulus_bitlen;
    }
  while (!gmpmee_done(t, test_time));

}

/* LCOV_EXCL_START */
void
usage(char *command_name) {
  printf("Usage: %s <ms>\n", command_name);
}
/* LCOV_EXCL_STOP */

int
main(int args, char *argv[])
{
  long ms;

  /* LCOV_EXCL_START */
  if (args == 2)
    {
      if (sscanf(argv[1], "%ld", &ms) != 1 || ms <= 0 || 60000 <= ms) {
        fprintf(stderr, "Not an integer! (%s)\n", argv[1]);
        exit(1);
      }
    }
  else
    {
      usage(argv[0]);
      exit(0);
    }
  /* LCOV_EXCL_STOP */

  printf("\n================================================================\n");
  printf("\n          TESTING GMPMEE\n\n");
  printf("================================================================\n\n");

  printf("Testing simultaneous exponentiation (%ld ms)... ", ms);
  test_spowm(ms);
  printf("done.\n");

  printf("Testing fixed base exponentiation (%ld ms)... ", ms);
  test_fpowm(ms);
  printf("done.\n");

  printf("Testing Miller-Rabin (%ld ms)... ", ms);
  test_miller_rabin(0, ms);
  printf("done.\n");

  printf("Testing Miller-Rabin safe prime (%ld ms)... ", ms);
  test_miller_rabin(1, ms);
  printf("done.\n");

  printf("Testing Miller-Rabin next (%ld ms)... ", ms);
  test_miller_rabin(2, ms);
  printf("done.\n");

  printf("Testing Miller-Rabin next safe prime (%ld ms)... ", ms);
  test_miller_rabin(3, ms);
  printf("done.\n\n");

  exit(0);
}
