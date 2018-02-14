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


/**
 *
 * Copyright 2008 2009 2010 2011 2013 Torbj&ouml;rn Granlund, Douglas
 * Wikstr&ouml;m
 *
 * <p>
 *
 * This file is part of GMP Modular Exponentiation Extension (GMPMEE).
 *
 * <p>
 *
 * GMPMEE is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * <p>
 *
 * GMPMEE is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * <p>
 *
 * You should have received a copy of the GNU General Public License
 * along with GMPMEE. If not, see <http://www.gnu.org/licenses/>.
 *
 * <hr>
 *
 * This is a minor extension of the <a
 * href="http://www.gmplib.org">Gnu Multiprecision Library (GNU
 * MP)</a>. It adds simultaneous modular exponentiation and fixed base
 * modular exponentiation functionality to the set of integer
 * functions.
 */



/**
 * @file gmpmee.h
 *
 * Implementation of simultaneous and fixed base modular
 * exponentiation using GMP routines, as well as (safe) primality
 * testing routines.
 *
 * <p>
 *
 * Simultaneous exponentiation is implemented in the usal way, see
 * <i>Handbook of Applied Cryptography</i>, Menezes, Oorshot, and
 * Vanstone, with a minor twist. The generators
 * <i>g<sub>1</sub>,...,g<sub>N</sub></i> are first divided into
 * <i>batches</i> of a user defined number of generators. Each such
 * batch is processed independently and the result is the product of
 * these partial results. Assume that all generators are in the same
 * batch. Then the generators are further subdivided into
 * <i>blocks</i>. For each block we precompute all possible products
 * of the generators of the block. When processing the batch, all
 * squarings are done for the complete batch, but multiplications are
 * done in the usual table driven way for each block.
 *
 * <p>
 *
 * Fixed base exponentiation is translated into simultaneous
 * exponentiation by considering the bases <i>g<sub>0</sub>=g</i>,
 * <i>g<sub>1</sub>=g<sup>2<sup>b</sup></sup></i>,
 * <i>g<sub>2</sub>=g<sup>2<sup>2b</sup></sup></i>,... Then when
 * computing <i>g<sup>x</sup></i>, the exponent is split into
 * <i>(x<sub>0</sub>,x<sub>1</sub>,...)</i>, where <i>0 <=
 * x<sub>i</sub> < 2<sup>b</sup></i>, and the result is given by:
 *
 * <p>
 *
 * <i>g<sub>1</sub><sup>x<sub>0</sub></sup>
 * g<sub>2</sub><sup>x<sub>2</sub></sup>
 * g<sub>3</sub><sup>x<sub>3</sub></sup>...</i>
 *
 * <p>
 *
 * The value of <i>b</i> is chosen for a given expected bit length of
 * <i>x</i>.
 *
 */

#ifndef GMPMEE_H
#define GMPMEE_H

#include <stdio.h>
#include <gmp.h>

/**
 * We use compiler flags that enforce that unused variables are
 * flagged as errors. In some cases where we are forced to keep
 * parameters to satisfy a given parameter signature, or when the
 * compiler complains incorrectly, this solves this in a documented
 * way.
 */
#define GMPMEE_UNUSED(x) ((void)(x))


/* #################### Simultaneous Exponentiation #################### */


/**
 * Stores the tables of precomputed products of subsets of the
 * bases. Each table contains the precomputed products for a range of
 * a given width of the bases.
 */
typedef struct
{
  size_t len;             /**< Total number of bases/exponents. */
  size_t block_width;     /**< Number of bases/exponents in each block. */
  size_t tabs_len;        /**< Number of blocks. */
  mpz_t **tabs;           /**< Table of tables, one sub-table for each block. */
  mpz_t modulus;          /**< Modulus used in computations. */

} gmpmee_spowm_tab[1]; /* Magic references. */


/**
 * Allocates and initializes a table for the given modulus, block
 * width, and total number of bases.
 *
 * @param table Table to be initialized
 * @param len Number of bases in the simultaneous exponentiation.
 * @param modulus Modulus.
 * @param block_width Number of bases used to build each subtable.
 */
void
gmpmee_spowm_init(gmpmee_spowm_tab table, size_t len, mpz_t modulus,
		  size_t block_width);

/**
 * Frees the memory allocated by table.
 *
 * @param table Table to be deallocated.
 */
void
gmpmee_spowm_clear(gmpmee_spowm_tab table);

/**
 * Fills the table with precomputed values using the given bases. The
 * array of bases must be of the length for which the table was
 * allocated.
 *
 * @param table Table to be initialized.
 * @param bases Bases for which precomputation is performed.
 */
void
gmpmee_spowm_precomp(gmpmee_spowm_tab table, mpz_t *bases);

/**
 * Computes a simultaneous exponentiation using the given table and
 * exponents. The number of exponents must match the number of bases
 * that was used during precomputation.
 *
 * @param rop Destination of result.
 * @param table Precomputed table representing the bases used.
 * @param exponents Exponents used in simultaneous exponentiation.
 */
void
gmpmee_spowm_table(mpz_t rop, gmpmee_spowm_tab table, mpz_t *exponents);

/**
 * Computes a simultaneous exponentiation. Precomputation is performed
 * in blocks of the given width in batches of the given batch size.
 *
 * @param rop Destination of result.
 * @param bases Bases for which precomputation is performed.
 * @param exponents Exponents used in simultaneous exponentiation.
 * @param len Number of bases in the simultaneous exponentiation.
 * @param modulus Modulus.
 * @param block_width Number of bases used to build each subtable.
 * @param batch_len Number of bases in each batch, where each batch
 * is computed independently.
 */
void
gmpmee_spowm_block_batch(mpz_t rop, mpz_t *bases, mpz_t *exponents, size_t len,
			 mpz_t modulus, size_t block_width, size_t batch_len);

/**
 * Computes a simultaneous exponentiation. Precomputation is performed
 * in blocks of a reasonable width in a single batch.
 *
 * @param rop Destination of result.
 * @param bases Bases for which precomputation is performed.
 * @param exponents Exponents used in simultaneous exponentiation.
 * @param len Number of bases in the simultaneous exponentiation.
 * @param modulus Modulus.
 */
void
gmpmee_spowm(mpz_t rop, mpz_t *bases, mpz_t *exponents, size_t len,
	     mpz_t modulus);

/**
 * Naively computes the exponentiated product of the bases to the
 * powers of the exponents modulo the given modulus. This is used for
 * debugging.
 *
 * @param rop Destination of result.
 * @param bases Bases.
 * @param exponents Exponents. used in simultaneous exponentiation.
 * @param len Number of bases.
 * @param modulus Modulus.
 */
void
gmpmee_spowm_naive(mpz_t rop, mpz_t *bases, mpz_t *exponents, size_t len,
		   mpz_t modulus);


/* #################### Fixed-Base Exponentiation #################### */

/**
 * Stores a fixed base exponentiation table.
 */
typedef struct
{
  gmpmee_spowm_tab spowm_table; /**< We exploit simultaneous exp. table. */
  size_t stretch;               /**< Normal number of bits of each
				   "subexponent". */
} gmpmee_fpowm_tab[1]; /* Magic references. */

/**
 * Allocates and initializes a table with the given modulus, block
 * width, and expected exponent bit length.
 *
 * @param table Table to be initialized
 * @param modulus Modulus.
 * @param block_width Number of bases used to build each subtable.
 * @param exponent_bitlen Expected bit length of exponent.
 */
void
gmpmee_fpowm_init(gmpmee_fpowm_tab table, mpz_t modulus,
		  size_t block_width, size_t exponent_bitlen);

/**
 * Frees the memory allocated by table.
 *
 * @param table Table to be deallocated.
 */
void
gmpmee_fpowm_clear(gmpmee_fpowm_tab table);

/**
 * Fills the table with precomputed values using the given basis.
 *
 * @param table Table to be initialized.
 * @param basis Basis for which precomputation is performed.
 */
void
gmpmee_fpowm_precomp(gmpmee_fpowm_tab table, mpz_t basis);

/**
 * Equivalent to calling gmpmee_fpowm_init and then gmpmee_fpowm_precomp.
 *
 * @param table Table to be initialized
 * @param basis Basis for which precomputation is performed.
 * @param modulus Modulus.
 * @param block_width Number of bases used to build each subtable.
 * @param exponent_bitlen Expected bit length of exponent.
 */
void
gmpmee_fpowm_init_precomp(gmpmee_fpowm_tab table, mpz_t basis,
			  mpz_t modulus, size_t block_width,
			  size_t exponent_bitlen);

/**
 * Computes a fixed base exponentiation using the given table and
 * exponent.
 *
 * @param rop Destination of result.
 * @param table Precomputed table representing the basis used.
 * @param exponent Exponent.
 */
void
gmpmee_fpowm(mpz_t rop, gmpmee_fpowm_tab table, mpz_t exponent);



/* #################### Primality Testing #################### */

/**
 * Stores state inbetween individual invokations of the Miller-Rabin
 * test and keeps allocated space.
 */
typedef struct
{
  mpz_t n;             /**< integer to be tested */
  mpz_t n_minus_1;     /**< n minus one */
  mpz_t q;             /**< q is defined by n=q*2^k+1 */
  unsigned long int k; /**< k is defined by n=q*2^k+1 */
  mpz_t y;             /**< y is temporary space */
} gmpmee_millerrabin_state[1]; /* Magic references. */


/**
 * Allocate and initialize Miller-Rabin state using the given integer.
 *
 * @param state State for testing.
 * @param n Integer to test.
 */
void
gmpmee_millerrabin_init(gmpmee_millerrabin_state state, mpz_t n);

/**
 * Updates the state to correspond to the next larger candidate
 * integer that passes the trial divisions.
 *
 * @param state State for testing primality.
 */
void
gmpmee_millerrabin_next_cand(gmpmee_millerrabin_state state);

/**
 * Free memory resources allocated for testing.
 *
 * @param state State for testing.
 */
void
gmpmee_millerrabin_clear(gmpmee_millerrabin_state state);

/**
 * Performs trial divisions and returns 0 or 1 depending on if a small
 * factor of the integer has been found or not. Assumes that the input
 * is greater than three.
 *
 * @param n Integer to test.
 */
int
gmpmee_millerrabin_trial(mpz_t n);

/**
 * Executes one round of the Miller-Rabin test and returns 0 or 1
 * depending on if the tested integer is deemed to be composite or
 * not.
 *
 * @param state State for testing.

 * @param base Base element used for testing. This must be an integer
 * in [2,n-2].
 */
int
gmpmee_millerrabin_once(gmpmee_millerrabin_state state, mpz_t base);

/**
 * Executes the Miller-Rabin test using randomness from one of GMP's
 * random sources. Assumes that the tested integer is greater than
 * three.
 *
 * @param rstate Source of randomness.
 * @param state State for testing.
 * @param reps Number of repetitions.
 */
int
gmpmee_millerrabin_reps_rs(gmp_randstate_t rstate,
			   gmpmee_millerrabin_state state,
			   int reps);

/**
 * Executes a number or repetitions of the Miller-Rabin test using
 * basis derived from the given GMP random source and returns 0 or 1
 * depending on if the tested integer is deemed to be composite or
 * not.
 *
 * <p>
 *
 * WARNING! GMP's random number generators are NOT cryptographically
 * secure.
 *
 * <p>
 *
 * @param rstate State of random number generator.
 * @param n Integer to test.
 * @param reps Repetitions of the Miller-Rabin test performed.
 */
int
gmpmee_millerrabin_rs(gmp_randstate_t rstate, mpz_t n, int reps);

/**
 * Searches for the smallest prime larger than the given
 * integer. Primality testing is done using the Miller-Rabin test
 * using randomness from one of GMP's random sources.
 *
 * @param rop Result destination.
 * @param rstate Source of randomness.
 * @param n Starting point in search.
 * @param reps Number of repetitions.
 */
void
gmpmee_millerrabin_next_rs(mpz_t rop, gmp_randstate_t rstate,
			   mpz_t n, int reps);

/**
 * Stores the states needed for using the Miller-Rabin test for
 * testing for safe-primality.
 */
typedef struct
{
  /**
   * State of the integer <i>n</i> to be tested.
   */
  gmpmee_millerrabin_state nstate;

  /**
   * State of the integer <i>(n-1)/2</i> to be tested.
   */
  gmpmee_millerrabin_state mstate;

} gmpmee_millerrabin_safe_state[1];

/**
 * Initialize Miller-Rabin state to be used for safe-primality
 * testing using the given integer.
 *
 * @param state State for testing safe-primality.
 * @param n Integer to test.
 */
void
gmpmee_millerrabin_safe_init(gmpmee_millerrabin_safe_state state, mpz_t n);

/**
 * Sets the state to the next candidate integer larger than the most
 * recently tested candidate that passes the trial divisions.
 *
 * @param state State for testing safe-primality.
 */
void
gmpmee_millerrabin_safe_next_cand(gmpmee_millerrabin_safe_state state);

/**
 * Free memory allocated in the states.
 *
 * @param state State for testing safe-primality.
 */
void
gmpmee_millerrabin_safe_clear(gmpmee_millerrabin_safe_state state);

/**
 * Performs trial divisions and returns 0 or 1 depending on if the
 * integer is definitely a not a safe prime, or if it could
 * potentially be a safe prime. Assumes that n is at least 8.
 *
 * @param n Integer to test.
 */

int
gmpmee_millerrabin_safe_trial(mpz_t n);

/**
 * Executes one round of the Miller-Rabin test and returns 0 or 1
 * depending on if the tested integer is deemed to not be a safe
 * prime, or a safe prime. Assumes that the tested integer is at least
 * 8.
 *
 * @param state State for testing safe-primality.
 * @param nbase Base element used for testing safe-primality. This
 * must be an integer in [2,n-1], where n is the integer to be tested.
 * @param mbase Base element used for testing safe-primality.  This
 * must be an integer in [2,m-1], where n=2m+1.
 */
int
gmpmee_millerrabin_safe_once(gmpmee_millerrabin_safe_state state,
			     mpz_t nbase, mpz_t mbase);

/**
 * Executes a safe primality test for the integer used to initialize
 * the given testing state, using randomness from the given GMP's
 * random source.
 *
 * @param rstate Source of randomness.
 * @param state State for testing safe-primality.
 * @param reps Number of repetitions.
 */
int
gmpmee_millerrabin_safe_reps_rs(gmp_randstate_t rstate,
				gmpmee_millerrabin_safe_state state,
				int reps);

/**
 * Executes several repetitions of the of the Miller-Rabin test and
 * returns 0 or 1 depending on if the tested integer is deemed to not
 * be a safe prime, or a safe prime. The basis elements are derived
 * from the given random number generator.
 *
 * <p>
 *
 * WARNING! GMP's random number generators are NOT cryptographically
 * secure.
 *
 * <p>
 *
 * @param rstate State of random number generator.
 * @param n Integer to test.
 * @param reps Repetitions of the Miller-Rabin test performed.
 */
int
gmpmee_millerrabin_safe_rs(gmp_randstate_t rstate, mpz_t n, int reps);

/**
 * Uses gmpmee_millerrabin_safe_rs to find the smallest safe prime
 * larger than the input integer.
 *
 * <p>
 *
 * WARNING! GMP's random number generators are NOT cryptographically
 * secure.
 *
 * <p>
 *
 * @param rop Found safe prime.
 * @param rstate State of random number generator.
 * @param n Integer to test.
 * @param reps Repetitions of the Miller-Rabin test performed.
 */
void
gmpmee_millerrabin_safe_next_rs(mpz_t rop, gmp_randstate_t rstate,
				mpz_t n, int reps);

/* #################### Utility Functions #################### */

/**
 * Naive implementation of a search for the next prime test.  Based on
 * GMP's probab_prime_p. Used for debugging.
 *
 * @param rop Next prime greater than the input integer.
 * @param n Starting point of search.
 * @param reps Number of repetitions of the Miller-Rabin test.
 */
void
mpz_probab_prime_p_next(mpz_t rop, mpz_t n, int reps);

/**
 * Naive implementation of a safe-primality test based on GMP's
 * probab_prime_p. Used for debugging.
 *
 * @param n Integer to be tested.
 * @param reps Number of repetitions of the Miller-Rabin test.
 */
int
mpz_probab_safe_prime_p(mpz_t n, int reps);

/**
 * Naive implementation of a search for the next safe prime test.
 * Based on probab_safe_prime_p. Used for debugging.
 *
 * @param rop Next safe prime greater than the input integer.
 * @param n Starting point of search.
 * @param reps Number of repetitions of the Miller-Rabin test.
 */
void
mpz_probab_safe_prime_p_next(mpz_t rop, mpz_t n, int reps);


/* #################### Utility Functions #################### */

/**
 * Allocates an array of <code>len</code> <code>mpz_t</code>.
 *
 * @param len Number of elements in array.
 * @return Pointer to allocated array.
 */
mpz_t*
gmpmee_array_alloc(size_t len);

/**
 * Allocates and initializes an array of <code>len</code>
 * <code>mpz_t</code>.
 *
 * @param len Number of elements in array.
 * @return Pointer to allocated array.
 */
mpz_t*
gmpmee_array_alloc_init(size_t len);

/**
 * Clears and deallocates the array containing <code>len</code>
 * <code>mpz_t</code>.
 *
 * @param a Array to be cleared and deallocated.
 * @param len Number of elements in array.
 */
void
gmpmee_array_clear_dealloc(mpz_t *a, size_t len);

/**
 * Fills the array rop containing <code>len</code> <code>mpz_t</code>
 * with random positive <code>n</code>-bit integers. <b>WARNING! The
 * pseudo-random generator of GMP used as a subroutine is *not*
 * cryptographically secure.</b>
 *
 * @param rop Destination of result.
 * @param len Number of elements in array.
 * @param state State of pseudo-random generator.
 * @param n Number of bits in each random integer.
 */
void
gmpmee_array_urandomb(mpz_t *rop, size_t len, gmp_randstate_t state,
		      unsigned long int n);

#endif /* GMPMEE_H */
