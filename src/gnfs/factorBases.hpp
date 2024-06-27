#ifndef FACTOR_BASES_HPP
#define FACTOR_BASES_HPP

#include "gnfs.hpp"     // serve perché ha la definizione di "ideale" (= first degree prime ideal o una potenza)

std::pair<uint8_t, bool> genPrimesList(const ZZ& n, primeList& primes, const ZZ& B); // se un primo divide n, .second è true e il divisore è in .back
std::pair<long, bool> buildQCB(const ZZ& n, factorBase& QCB, const ZZX&f, const ZZ& L, long t);
void rootsOfFmodP(const ZZX& f, const ZZ p, std::vector<ZZ>& roots, bool multiple=true);
long buildFBase(factorBase& FB, const ZZX&f, const ZZ& B, primeList& primes);

#endif