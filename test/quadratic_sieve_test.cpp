#include <gtest/gtest.h>

#include "../src/quadratic_sieve/quadratic_sieve.hpp"

#include <gmp.h>
#include <gmpxx.h>
#include <tuple>
#include <vector>

#include "../src/utils/utils_mpz/hash_mpz.hpp"
#include "../src/utils/utils_mpz/mpz_ull.hpp"
#include "../src/utils/utils_modP/roots_modP.hpp"


TEST(Quadratic, NeqPQ) {
  mpz_class n = 700757277917;
  mpz_class p, q;
  quadratic_sieve(n, p, q);

  ASSERT_EQ(n, p*q);
}



// controllare che i vari componenti funzionino

TEST(Quadratic, factorBase){
  mpz_class n, B;
  n = "5137851827";
  B = 383;
  std::vector<std::pair<mpz_class, unsigned short>> factorBase;
  buildFactorBase(n, B, factorBase);

  ASSERT_GT(factorBase.size(), 2);

  std::vector<unsigned long long> primi1, primi2 = {2, 11, 13, 17, 19, 23, 37, 41, 43, 59, 73, 103, 107, 113, 127, 131, 137, 149, 163, 179, 181,
191, 197, 199, 211, 223, 233, 241, 251, 271, 281, 293, 307, 311, 317, 347, 367, 373, 379, 383};

  mpz_class P, tmp1; unsigned short l;
  for(auto coppia : factorBase){
    std::tie(P, l) = coppia;

    int probPrime = mpz_probab_prime_p(P.get_mpz_t(), 20);
    EXPECT_GT(probPrime, 0) << "element " << P << " in factor base is not prime";

    double log2pD = log2(P.get_d());
    unsigned short l2 = (unsigned short) (log2pD * (1<<6));
    EXPECT_EQ(l, l2) << "logarithm in base 2 of " << P <<" is wrong";

    unsigned long long p = mpz_2_ull(P, tmp1);

    primi1.push_back(p);
  }

  //ASSERT_THAT(primi1, ContainerEq(primi2));
  ASSERT_EQ(primi1, primi2);
}

TEST(Quadratic, findRootModP3mod4){
  mpz_class root = 302, p = 7883, a = (root*root)%p, tmp1, tmp2, tmp3, v;
  findRoot(a, p, v, tmp1, tmp2, tmp3);
  ASSERT_EQ((v*v)%p, a);
}


TEST(Quadratic, findRootModP1mod4){
  mpz_class root = 102, p = 257, a = (root*root)%p, tmp1, tmp2, tmp3, v;
  findRoot(a, p, v, tmp1, tmp2, tmp3);
  ASSERT_EQ((v*v)%p, a);
}


TEST(Quadratic, initializedSive){
  mpz_class n, L;
  n = "758216904857";
  L = "20871";

  mpz_class base = sqrt(n) + 1;

  std::vector<std::pair<mpz_class, unsigned short>> factorBase = {
    {2, 64},    {11, 221},  {13, 236},  {19, 271},  {23, 289},  {29, 310},  {37, 333},  {53, 366},  {71, 393},  {73, 396},
    {79, 403},  {83, 408},  {89, 414},  {103, 427}, {109, 433}, {113, 436}, {127, 447}, {131, 450}, {137, 454}, {157, 466},
    {167, 472}, {179, 478}, {181, 479}, {193, 485}, {197, 487}, {199, 488}, {211, 494}, {223, 499}, {271, 517}, {281, 520},
    {283, 521}, {313, 530}, {317, 531}
  };

  std::unordered_map<mpz_class, elemSetaccio> setaccio;
  initializeSieve(n, base, L, factorBase, setaccio);

  // verificare che il setaccio sia inizializzato correttamente:
  // v- ogni fattore deve comparire in esattamente 2 punti (eccetto per 2 (1 punto), e le potenze di 2 da 8 in su (4 punti))
  // v- se un fattore P compare nell'entry di a, allora P|y(x)=x^2-n=(a+base)^2-n
  // - se P compare nelle entry a e b, allora P!|(a-b)
  // - non ho voglia di controllare i log

  std::unordered_map<mpz_class, std::vector<mpz_class>> occorrenze;
  mpz_class a, P, baseSquaredMinusN = base*base - n, twoBase = base<<1;   // y(x) = y(a+base) = baseSquaredMinusN + a*(twoBase+a)

  for(auto it = setaccio.begin(); it != setaccio.end(); it++){
    a = it->first;

    for(auto j = it->second.begin(); j != it->second.end(); j++){
      P = j->first;
      EXPECT_EQ((baseSquaredMinusN + a*(twoBase+a))%P, 0) << P << "non divide y(" << a+base << ") = " << (baseSquaredMinusN + a*(twoBase+a))%P;

      auto k = occorrenze.find(P);
      if(k == occorrenze.end()) occorrenze[P] = std::vector<mpz_class>{a};
      else occorrenze[P].push_back(a);
    }
  }

  for(auto it = occorrenze.begin(); it != occorrenze.end(); it++){
    P = it->first;

    if(P==2){
      ASSERT_GT(it->second.size(), 0) << "non ci sono radici di n=" << n << "modulo 2" ;
      ASSERT_LT(it->second.size(), 2) << "ci sono troppe radici di n=" << n << "modulo 2" ;
    } else if(mpz_even_p(P.get_mpz_t()) && P >= 8){
      ASSERT_GT(it->second.size(), 3) << "non ci sono abbastanza radici di n=" << n << "modulo " << P ;
      ASSERT_LT(it->second.size(), 5) << "ci sono troppe radici di n=" << n << "modulo " << P ;

      ASSERT_NE(((it->second)[0] - (it->second)[1])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
      ASSERT_NE(((it->second)[0] - (it->second)[2])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
      ASSERT_NE(((it->second)[0] - (it->second)[3])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
      ASSERT_NE(((it->second)[1] - (it->second)[2])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
      ASSERT_NE(((it->second)[1] - (it->second)[3])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
      ASSERT_NE(((it->second)[2] - (it->second)[3])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
    } else{
      ASSERT_GT(it->second.size(), 1) << "non ci sono abbastanza radici di n=" << n << "modulo " << P ;
      ASSERT_LT(it->second.size(), 3) << "ci sono troppe radici di n=" << n << "modulo " << P ;

      ASSERT_NE(((it->second)[0] - (it->second)[1])%P, 0) << "ci sono radici uguali di n=" << n << "modulo " << P;
    }
  }

}



/*TEST(Quadratic, bigInput) {
  mpz_class n;
  n = "2694510496740314556501037319858087";
  mpz_class p, q;
  quadratic_sieve(n, p, q);

  ASSERT_EQ(n, p*q);
}*/