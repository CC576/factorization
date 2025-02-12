#include <gtest/gtest.h>

#include "../src/trial_div/trial_division.hpp"

#include <gmp.h>
#include <gmpxx.h>
#include <tuple>

TEST(TrialDivision, NeqPQ) {
  mpz_class n = 700757277917;
  mpz_class p, q;
  trial_division(n, p, q);

  ASSERT_EQ(n, p*q);
}
