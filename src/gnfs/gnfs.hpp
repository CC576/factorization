#ifndef GENERALNUMBERFIELDSIEVE_HPP
#define GENERALNUMBERFIELDSIEVE_HPP

#include <gmp.h>
#include <gmpxx.h>

#include <blanczos.h>

#include <cmath>
#include <utility>
#include <vector>
#include <unordered_map>
//#include <map>
#include <bitset>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/pair_ZZX_long.h>
#include <NTL/ZZXFactoring.h>
//#include <NTL/ZZ_pE.h>
#include "../utils/utils_NTL/hash_ZZ.hpp"
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
struct smoothElemGNFS{
    ZZ a; ZZ b;
    std::vector<uint32_t> rationalPpos;     // indici dei primi razionali che dividono a-bm; include -1
    std::vector<uint32_t> algebraicPpos;    // indici dei first degree prime ideals che dividono <a-b*theta>; include -1 (come divisore di N(a-b*theta))
                                            // gli indici sono uint32_t per blanczos

    smoothElemGNFS(ZZ _a, ZZ _b){
        a = _a; b = _b;
    }
};

typedef std::unordered_map<ZZ, long> factorization;    // primo-esponente; non conta -1
typedef std::vector<std::pair<long, uint8_t>> primeList;
typedef std::vector<ideal> factorBase;
typedef std::vector<uint8_t> sieveArray;
typedef std::unordered_map<std::pair<ZZ, ZZ>, uint32_t> idealMap;

void gnfs(const mpz_class&, mpz_class&, mpz_class&);

void chooseParams(const mpz_class& n, long& d, ZZ& m, ZZX& f, ZZ& B);
bool findEarlyFactors(const ZZ& n, ZZ& fattore, const ZZX& f, ZZX& fPrime, const ZZ& m);
std::pair<uint8_t, bool> buildFactorBases(const ZZ& n, const ZZX& f, factorBase& RFB, factorBase& AFB, factorBase& QCB, const ZZ& B, ZZ& L, ZZ& m, primeList& primes, long& k, long& l, long& t);   // se .second è true, L contiene un fattore primo di N
std::pair<uint64_t, ZZ> sieve(const ZZ& n, const ZZ& m, const ZZX& f, const uint8_t logMaxP2, const ZZ& L, ZZ& b, const long maxA, const uint8_t logm, const long quadChars, long& numSmooths, long& rationalPrimes, long& algebraicPrimes,  const primeList& primes, const factorBase& RFB, const factorBase& AFB, std::vector<smoothElemGNFS>& smooths);   // restituisce il numero di entries
uint64_t getMatrix(std::vector<uint32_t>& mat, std::vector<smoothElemGNFS>& smooths, long numRatPrimes, long numAlgPrimes, const factorBase& QCB);

#endif