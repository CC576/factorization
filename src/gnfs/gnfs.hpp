#ifndef GENERALNUMBERFIELDSIEVE_HPP
#define GENERALNUMBERFIELDSIEVE_HPP

#include <gmp.h>
#include <gmpxx.h>

#include <blanczos.h>

#include <cmath>
#include <utility>
#include <vector>
//#include <bitset>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/pair_ZZX_long.h>
#include <NTL/ZZXFactoring.h>
using namespace NTL;

struct ideal{                           // first degree prime ideal o una potenza
    ZZ p;                               // p rappresenta "l'incremento": può essere p primo anche se si tratta di una potenza di ideale (in quel caso logP = e*log_2(p))
    ZZ r;
    uint8_t logP;

    ideal(const ZZ& _p, const ZZ& _r, const uint8_t _logP){
        p = _p; r = _r;
        logP = _logP;
    }
};

typedef std::vector<ideal> factorBase;

void gnfs(const mpz_class&, mpz_class&, mpz_class&);

void chooseParams(const mpz_class& n, long& d, ZZ& m, ZZX& f, ZZ& B);
bool findEarlyFactors(const ZZ& n, ZZ& fattore, const ZZX& f, ZZX fPrime, const ZZ& m);
uint8_t buildFactorBases(const ZZ& n, const ZZX& f, factorBase& RFB, factorBase& AFB, factorBase& QCB, const ZZ& B, ZZ& L, ZZ& m, std::vector<std::pair<long, uint8_t>>& primes);

#endif