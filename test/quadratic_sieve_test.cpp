#include <gtest/gtest.h>

#include "../src/quadratic_sieve/quadratic_sieve.hpp"

#include <gmp.h>
#include <gmpxx.h>
#include <tuple>
#include <vector>

#include "../src/utils/utils_mpz/hash_mpz.hpp"
#include "../src/utils/utils_mpz/mpz_ull.hpp"


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





/*TEST(Quadratic, bigInput) {
  mpz_class n;
  n = "2694510496740314556501037319858087";
  mpz_class p, q;
  quadratic_sieve(n, p, q);

  ASSERT_EQ(n, p*q);
}*/