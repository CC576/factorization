#ifndef MPZ_LL_HPP
#define MPZ_LL_HPP

#include <gmp.h>
#include <gmpxx.h>


unsigned long long mpz_2_ull(mpz_t&, mpz_t&);
unsigned long long mpz_2_ull(mpz_class&, mpz_class&);

void mpz_from_ull(unsigned long long, mpz_t&);
void mpz_from_ull(unsigned long long, mpz_class&);

#endif