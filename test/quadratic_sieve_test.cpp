#include <gtest/gtest.h>

#include "../src/quadratic_sieve/quadratic_sieve.hpp"

#include <gmp.h>
#include <gmpxx.h>
#include <tuple>

TEST(Fermat, NeqPQ) {
  mpz_class n = 700757277917;
  mpz_class p, q;
  quadratic_sieve(n, p, q);

  ASSERT_EQ(n, p*q);
}


TEST(Fermat, bigInput) {
  mpz_class n;
  n = "17070404940041398912663602465930233706123576130475774484108750122669";
  mpz_class p, q;
  quadratic_sieve(n, p, q);

  ASSERT_EQ(n, p*q);
}


//to-do: controllare che i vari componenti funzionino