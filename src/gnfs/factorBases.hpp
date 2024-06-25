#ifndef FACTOR_BASES_HPP
#define FACTOR_BASES_HPP

#include "gnfs.hpp"     // serve perché ha la definizione di "ideale" (= first degree prime ideal o una potenza)

uint8_t genPrimesList(std::vector<std::pair<long, uint8_t>>& primes, const ZZ& B);
long buildQCB(factorBase& QCB, const ZZX&f, const ZZ& L, long t);
void rootsOfFmodP(const ZZX& f, const ZZ p, std::vector<ZZ>& roots, bool multiple=true);
long buildFBase(factorBase& FB, const ZZX&f, const ZZ& B, std::vector<std::pair<long, uint8_t>>& primes);

#endif