#include <gtest/gtest.h>

#include "../src/fermat/fermat.hpp"

#include <gmp.h>
#include <gmpxx.h>
#include <tuple>

TEST(Fermat, NeqPQ) {
  mpz_class n = 700757277917;
  mpz_class p, q;
  fermat(n, p, q);

  ASSERT_EQ(n, p*q);
}