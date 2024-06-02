#ifndef ROOTS_MODP_HPP
#define ROOTS_MODP_HPP

#include <gmp.h>
#include <gmpxx.h>

#include "../utils_mpz/mpz_ull.hpp"

void findRoot(const mpz_class&, const mpz_class&, mpz_class&, mpz_class&, mpz_class&, mpz_class&);
void getNonQResidue(const mpz_class&, mpz_class&);
unsigned long getOrdInSylow(const mpz_class&, mpz_class&, unsigned long, mpz_class&);
void extendRoot2(mpz_class&, mpz_class&, mpz_class&, mpz_class&);
void extendRootP(mpz_class&, mpz_class&, mpz_class&, mpz_class&, mpz_class&, mpz_class&, mpz_class&);

#endif