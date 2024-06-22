#include <gtest/gtest.h>

#include "../src/gnfs/gnfs.hpp"

#include <gmp.h>
#include <gmpxx.h>

//#include <NTL/tools.h>
//using namespace NTL;

//#include <tuple>
//#include <vector>
//#include <algorithm>

//#include "../src/utils/utils_mpz/hash_mpz.hpp"
//#include "../src/utils/utils_mpz/mpz_ull.hpp"
//#include "../src/utils/utils_modP/roots_modP.hpp"


TEST(Quadratic, NeqPQ) {
  mpz_class n = 700757277917;
  mpz_class p, q;
  gnfs(n, p, q);

  ASSERT_EQ(n, p*q);
}



// controllare che i vari componenti funzionino














// test per blanczos? O assumo che funzioni?


TEST(GNFS, bigInput) {
  mpz_class n;
  n = "2694510496740314556501037319858087";   // 100710366161784104467949462748056679263
  mpz_class p, q;
  gnfs(n, p, q);

  ASSERT_EQ(n, p*q);
}