#ifndef MPZ_ZZ_HPP
#define MPZ_ZZ_HPP

#include <gmp.h>
#include <gmpxx.h>
#include <NTL/ZZ.h>
using namespace NTL;

void mpz_2_ZZ(const mpz_t&, ZZ&);
void mpz_2_ZZ(const mpz_class&, ZZ&);

void mpz_from_ZZ(const ZZ&, mpz_t&);
void mpz_from_ZZ(const ZZ&, mpz_class&);

#endif