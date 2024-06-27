#include "gnfs.hpp"

void computeFactorization(const ZZ& m, const ZZX& f, const primeList& primes, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, factorization& fattorizzazione, bool algebraic=false);
void computeProdMod(const ZZ& modulo, const factorization& fat, ZZ_p& res);
void estimateXSize(const ZZX& f, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, ZZ& res);
void findNextInertPrime(const ZZX& f, const ZZ& last, ZZ& inertPrime);
void compSquareRootInGF(const ZZX& f, const ZZ& l, const ZZX& fPrime, const factorization& algFat, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, ZZ_pX& res);
