#ifndef GENERALNUMBERFIELDSIEVE_HPP
#define GENERALNUMBERFIELDSIEVE_HPP

#include <gmp.h>
#include <gmpxx.h>

#include <cmath>
#include <blanczos.h>
//#include <bitset>

#include <NTL/tools.h>
using namespace NTL;

void gnfs(mpz_class&, mpz_class&, mpz_class&);

#endif