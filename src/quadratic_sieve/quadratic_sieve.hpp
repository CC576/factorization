#ifndef QUADRATIC_HPP
#define QUADRATIC_HPP

#include <gmp.h>
#include <gmpxx.h>
//#include <mpfr.h>
#include <cmath>
#include <blanczos.h>
#include <utility>
#include <tuple>

void quadratic_sieve(mpz_class&, mpz_class&, mpz_class&);

void choose_params(mpz_class&, mpz_class&, mpz_class&);

#endif