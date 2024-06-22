#ifndef GENERALNUMBERFIELDSIEVE_HPP
#define GENERALNUMBERFIELDSIEVE_HPP

#include <gmp.h>
#include <gmpxx.h>

#include <cmath>
#include <blanczos.h>
//#include <bitset>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
using namespace NTL;

void gnfs(const mpz_class&, mpz_class&, mpz_class&);

void chooseParams(const mpz_class& n, long& d, ZZ& m, ZZX& f, ZZ& B);

#endif