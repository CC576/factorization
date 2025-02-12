#ifndef SIEVING_HPP
#define SIEVING_HPP

#include "gnfs.hpp"

void initSeed();
long lineSieve(const ZZX& f, const uint8_t logMaxP2, const ZZ& b, const long maxA, const uint8_t logm, const factorBase& RFB, const factorBase& AFB, sieveArray& setaccio, std::vector<bool>& isProbSmooth, std::vector<ZZ>& probSmooths);
std::pair<bool, bool> trialDivide(const ZZ& a, const ZZ &b, const ZZ& num, const ZZ& largePrimeBound, const primeList& primes, std::vector<std::pair<ZZ, ZZ>>& fattori);     // dice se l'elemento era smooth (anche partialmente) o no e totalmente smooth o parziale
void moveFactors(std::vector<uint32_t>& indici, std::vector<std::pair<ZZ, ZZ>>& fattori, const ZZ& maxP, ZZ& maxNewP, idealMap& usedPrimes, long& newPrimes);


uint64_t sieveForTesting(const ZZ& n, const ZZ& m, const ZZX& f, const uint8_t logMaxP2, const ZZ& L, ZZ& b, const long maxA, const uint8_t logm, const long quadChars, long& numSmooths, long& rationalPrimes, long& algebraicPrimes,  const primeList& primes, const factorBase& RFB, const factorBase& AFB, std::vector<smoothElemGNFS>& smooths, std::vector<ZZ>& pos2ratPrimes, std::vector<ZZ>& pos2algPrimes);

#endif