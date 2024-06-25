#include <gtest/gtest.h>

#include "../src/gnfs/gnfs.hpp"
#include "../src/gnfs/factorBases.hpp"

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


TEST(GNFS, rootsModP){
  ZZX f;
  SetCoeff(f, 0, 8);
  SetCoeff(f, 1, 29);
  SetCoeff(f, 2, 15);
  SetCoeff(f, 3,1);

  ZZ p(67);
  ZZ_p::init(p);
  //ZZ_pX fp = conv<ZZ_pX>(f);

  std::vector<ZZ> roots, trueRoots = {ZZ(2), ZZ(6), ZZ(44)};
  rootsOfFmodP(f, p, roots);

  sort(roots.begin(), roots.end());
  ASSERT_EQ(roots, trueRoots);
}



TEST(GNFS, factorBases){
  ZZ n(2601699059);       // questo n genera un polinomio f con radici multiple modulo almeno 4 primi p, di cui uno dà esponente 4
  mpz_class nMPZ = 2601699059;

  long d;
  ZZ m, B;
  ZZX f;
  chooseParams(nMPZ, d, m, f, B);

  factorBase RFB, AFB, QCB;
  ZZ L, tmp;
  std::vector<std::pair<long, uint8_t>> primes;
  uint8_t logMaxP2 = buildFactorBases(n, f, RFB, AFB, QCB, B, L, m, primes);

  ZZ p(2);//, last(2);
  //ZZ p, Pot, r;
  //uint8_t l, trueLog;

  ZZ_pX fp, q;
  ZZ_p tmpP;


  // check RFB
  for(auto& [Pot, r, l] : RFB){
    long e = 1;
    if(ProbPrime(Pot)){
      p = Pot;
    } else{ // assumo p = last, Pot potenza di p, per l'ordine in cui ho aggiunto gli elementi
      tmp = p;
      while(tmp < Pot){
        tmp*=p;
        e++;
      }
      ASSERT_EQ(tmp, Pot) << "p=" << Pot << "in RFB is neither prime nor a power of a prime";
    }
    // Pot = p^e

    uint8_t trueLog = (uint8_t) log2(conv<long>(p));  // per RFB e AFB p<B dovrebbe stare in un long
    if(trueLog == l){     // caso radice singola
      EXPECT_TRUE(m%Pot == r) << "m=" << m << " is not r=" << r << "mod Pot=" << Pot << " in RFB" ;
    } else{               // caso radice multipla
      ASSERT_TRUE(false) << "There shouldn't be multiple roots in RFB, how did we even get here?";
    }

  }


  // check AFB
  for(auto& [Pot, r, l] : AFB){
      long e = 1;
      if(ProbPrime(Pot)){
        p = Pot;
      } else{ // assumo p = last, Pot potenza di p, per l'ordine in cui ho aggiunto gli elementi
        tmp = p;
        while(tmp < Pot){
          tmp*=p;
          e++;
        }
        ASSERT_EQ(tmp, Pot) << "p=" << Pot << "in AFB is neither prime nor a power of a prime";
      }
      // Pot = p^e

      ZZ_p::init(Pot);
      fp = conv<ZZ_pX>(f);
      uint8_t trueLog = (uint8_t) log2(conv<long>(p));  // per RFB e AFB p<B dovrebbe stare in un long

      if(trueLog == l){     // caso radice singola

        tmpP = eval(fp, conv<ZZ_p>(r));
        EXPECT_TRUE(IsZero(tmpP)) << "r=" << r << " in AFB is not a root of f mod Pot=" << Pot;

      } else{               // caso radice multipla

        // per come ho calcolato l deve valere questa cosa
        EXPECT_TRUE(l % trueLog == 0) << "Log_2 l=" << l << " of power " << Pot << " = " << p <<"^"<<e<<" in AFB is wrong";
        e = l/trueLog;
        Pot = power(p, (long) e);

        for(tmp = 0; tmp < Pot; tmp+=p){
          tmpP = eval(fp, conv<ZZ_p>(r+tmp));
          ASSERT_TRUE(IsZero(tmpP)) << "r=" << r+tmp << " in AFB is not a root of f mod Pot=" << Pot;
        }
      }

    }


  // check QCB
  long t = log2(conv<double>(n));
  EXPECT_GE(QCB.size(), 3*t) << "Size QCB is too small";
  EXPECT_LE(QCB.size(), 3*(t+1)) << "Size QCB is too big";

  for(auto& [p, r, l] : QCB){
    ZZ_p::init(p);
    fp = conv<ZZ_pX>(f);

    ASSERT_TRUE(ProbPrime(p)) << "p=" << p << " in QCB is not prime";

    uint8_t trueLog = (uint8_t) log2(conv<double>(p));
    if(trueLog >= l) EXPECT_LT(trueLog - l, 1) << "Log2 of p=" << p <<" in QCB is wrong";
    else EXPECT_LT(l - trueLog, 1) << "Log2 of p=" << p <<" in QCB is wrong";

    SetX(q); q-=conv<ZZ_p>(r);      // q = x-r
    EXPECT_TRUE(divide(fp, q)) << "r=" << r << " in QCB is not a root of f mod p=" << p;
    q = q*q;
    EXPECT_FALSE(divide(fp, q)) << "r=" << r << " in QCB is a multiple root of f mod p=" << p;
  }

} // il test è stato superato, il mio calcolo sulle radici mod potenze di primi era (molto probabilmente) giusto!





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