#include <gtest/gtest.h>

#include "../src/gnfs/gnfs.hpp"

#include <gmp.h>
#include <gmpxx.h>

//#include <NTL/ZZ.h>
//using namespace NTL;

//#include <tuple>
//#include <vector>
//#include <algorithm>

//#include "../src/utils/utils_mpz/hash_mpz.hpp"
//#include "../src/utils/utils_mpz/mpz_ull.hpp"
//#include "../src/utils/utils_modP/roots_modP.hpp"
#include "../src/utils/utils_mpz/mpz_ZZ.hpp"
#include "../src/utils/utils_NTL/utils_ZZX.hpp"


TEST(GNFS, NeqPQ) {
  mpz_class n = 700757277917;
  mpz_class p, q;
  gnfs(n, p, q);

  ASSERT_EQ(n, p*q);
}



// controllare che i vari componenti funzionino

TEST(GNFS, fmEqN){
  mpz_class nMPZ = 700757277917;
  long d;
  ZZ m, B;
  ZZX f;
  chooseParams(nMPZ, d, m, f, B);
  ZZ n, fm;
  mpz_2_ZZ(nMPZ, n);
  ZZX_eval(fm, f, m);
  ASSERT_EQ(n, fm);
}


TEST(GNFS, EarlyFactors){
  ZZ p, q, n, m, tmp;
  p = 137; q = 157;
  n = p*q; m = 26;
  ZZX f1, f2, f3, f, tmpF;
  SetX(f1); SetX(f2); SetX(f3);
  SetCoeff(f1, 0, 1-m);
  SetCoeff(f2, 0, p-m);
  SetCoeff(f3, 0, q-m);
  f = f1*f2*f3;

  ZZX_eval(tmp, f, m);
  ASSERT_EQ(n, tmp);

  findEarlyFactors(n, tmp, f, tmpF, m);
  tmp = abs(tmp);
  ASSERT_TRUE(tmp == p || tmp == q);
}








// test per blanczos? O assumo che funzioni?


#ifndef DEBUG
TEST(GNFS, bigInput) {
  mpz_class n;
  n = "2694510496740314556501037319858087";   // 100710366161784104467949462748056679263
  mpz_class p, q;
  gnfs(n, p, q);

  ASSERT_EQ(n, p*q);
}
#endif