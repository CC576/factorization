#ifndef GENERALNUMBERFIELDSIEVE_HPP
#define GENERALNUMBERFIELDSIEVE_HPP

#include <gmp.h>
#include <gmpxx.h>

#include <cmath>
#include <blanczos.h>
//#include <bitset>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/pair_ZZX_long.h>
#include <NTL/ZZXFactoring.h>
using namespace NTL;

void gnfs(const mpz_class&, mpz_class&, mpz_class&);

void chooseParams(const mpz_class& n, long& d, ZZ& m, ZZX& f, ZZ& B);
bool findEarlyFactors(const ZZ& n, ZZ& fattore, const ZZX& f, ZZX fPrime, const ZZ& m);

#endif