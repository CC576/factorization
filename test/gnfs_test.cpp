#include <gtest/gtest.h>

#include "../src/gnfs/gnfs.hpp"
#include "../src/gnfs/factorBases.hpp"
#include "../src/gnfs/sieving.hpp"
#include "../src/gnfs/computeSquareRoots.hpp"

#include <gmp.h>
#include <gmpxx.h>

//#include <NTL/ZZ.h>
//using namespace NTL;
#include <NTL/ZZ_pXFactoring.h>

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

  long d, k, l, t;
  ZZ m, B;
  ZZX f;
  chooseParams(nMPZ, d, m, f, B);

  factorBase RFB, AFB, QCB;
  ZZ L, tmp;
  std::vector<std::pair<long, uint8_t>> primes;
  auto [logMaxP2, _] = buildFactorBases(n, f, RFB, AFB, QCB, B, L, m, primes, k, l, t);

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
  long t1 = log2(conv<double>(n));
  EXPECT_GE(QCB.size(), 3*t1) << "Size QCB is too small";
  EXPECT_LE(QCB.size(), 3*(t1+1)) << "Size QCB is too big";

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


TEST(GNFS, sieveResult){
  ZZ n(2601699059);       // questo n genera un polinomio f con radici multiple modulo almeno 4 primi p, di cui uno dà esponente 4
  mpz_class nMPZ = 2601699059;

  long d, k, l, t;
  ZZ m, B;
  ZZX f;
  chooseParams(nMPZ, d, m, f, B);

  factorBase RFB, AFB, QCB;
  ZZ L, tmp;
  std::vector<std::pair<long, uint8_t>> primes;
  auto [logMaxP2, _] = buildFactorBases(n, f, RFB, AFB, QCB, B, L, m, primes, k, l, t);


  long maxA;
  if(B > m/2) maxA = conv<long>(m/2);
  else maxA = conv<long>(B);
  uint8_t logm = log2(conv<double>(m));
  ZZ b(1);

  std::vector<ZZ> pos2ratPrimes, pos2algPrimes;
  std::vector<smoothElemGNFS> smooths;
  long numSmooths = 0, rationalPrimes = 0, algebraicPrimes = 0;
  uint64_t entries = 0ull;
  entries = sieveForTesting(n, m, f, logMaxP2, L, b, maxA, logm, t, numSmooths, rationalPrimes, algebraicPrimes, primes, RFB, AFB, smooths, pos2ratPrimes, pos2algPrimes);

  ZZ aMinusBm, normAB, P, tmp2;
  for(const auto& coppia : smooths){
    aMinusBm = coppia.a - coppia.b*m;
    normAB = 0; tmp = 1;
    for(long i=0; i<=deg(f); i++){
      //std::cout << i << " " << coeff(f, i) << " " << tmp << " " << power(b, deg(f)-i) << std::endl;
      normAB += coeff(f, i)*tmp*power(coppia.b, deg(f)-i);
      tmp*=coppia.a;      // lo moltiplico con un accumulatore perché potrebbe essere 0
      //std::cout << normAB << " " << std::endl;
    }

    // assicurarsi che sia smooth (almeno parziale)
      // a-bm
      tmp = aMinusBm;
      ASSERT_FALSE(IsZero(tmp)) << "a-bm=" << tmp << "for (" << coppia.a << ", " << coppia.b << ")";
      tmp = abs(tmp);
      for(auto& primo : primes){
        P = conv<ZZ>(primo.first);
        while(tmp%P == 0) tmp/=P;
      }
      ASSERT_LE(tmp, L) << "(" << coppia.a << ", " << coppia.b << ") is not rational smooth";

      // N(a, b)
      tmp = normAB;
      ASSERT_FALSE(IsZero(tmp)) << "N(a,b)=" << tmp << "for (" << coppia.a << ", " << coppia.b << ")";
      tmp = abs(tmp);
      for(auto& primo : primes){
        P = conv<ZZ>(primo.first);
        while(tmp%P == 0) tmp/=P;
      }
      ASSERT_LE(tmp, L) << "(" << coppia.a << ", " << coppia.b << ") is not algebraic smooth";


    // sia per a-bm che per N(a,b) assicurarsi che tolti i fattori in rationalPpos e algebraicPpos restino dei quadrati
      // a-bm
      tmp = aMinusBm;
      for(auto& pos : coppia.rationalPpos){
        ASSERT_LT(pos, pos2ratPrimes.size()) << "Rational prime out of range";
        P = pos2ratPrimes[pos];
        if(P == -1) ASSERT_LT(tmp, 0) << "a-bm=" << aMinusBm << " is not negative for (" << coppia.a << ", " << coppia.b << ")";
        else ASSERT_TRUE(tmp%P == 0) << "a-bm=" << aMinusBm << " is not divisible by one of its factors for (" << coppia.a << ", " << coppia.b << ")";
        tmp /= P;
      }
      tmp2 = SqrRoot(tmp);
      ASSERT_EQ(tmp, tmp2*tmp2) << "The odd exponent factors of a-bm=" << aMinusBm << " are wrong for (" << coppia.a << ", " << coppia.b << ")";

      // N(a, b)
      tmp = normAB;
      for(auto& pos : coppia.algebraicPpos){
        ASSERT_LT(pos, pos2algPrimes.size()) << "Algebraic prime out of range";
        P = pos2algPrimes[pos];
        if(P == -1) ASSERT_LT(tmp, 0) << "N(a, b)=" << normAB << " is not negative for (" << coppia.a << ", " << coppia.b << ")";
        else ASSERT_TRUE(tmp%P == 0) << "N(a, b)=" << normAB << " is not divisible by one of its factors for (" << coppia.a << ", " << coppia.b << ")";
        tmp /= P;
      }
      tmp2 = SqrRoot(tmp);
      ASSERT_EQ(tmp, tmp2*tmp2) << "The odd exponent factors of N(a, b)=" << normAB << " are wrong for (" << coppia.a << ", " << coppia.b << ")";
  }

}





// test per blanczos? O assumo che funzioni?
TEST(GNFS, linearDependencies){
  ZZ n(2601699059);       // questo n genera un polinomio f con radici multiple modulo almeno 4 primi p, di cui uno dà esponente 4
  mpz_class nMPZ = 2601699059;

  long d, k, l, t;
  ZZ m, B;
  ZZX f;
  chooseParams(nMPZ, d, m, f, B);

  factorBase RFB, AFB, QCB;
  ZZ L, tmp;
  std::vector<std::pair<long, uint8_t>> primes;
  auto [logMaxP2, _] = buildFactorBases(n, f, RFB, AFB, QCB, B, L, m, primes, k, l, t);


  long maxA;
  if(B > m/2) maxA = conv<long>(m/2);
  else maxA = conv<long>(B);
  uint8_t logm = log2(conv<double>(m));
  ZZ b(1);

  std::vector<ZZ> pos2ratPrimes, pos2algPrimes;
  std::vector<smoothElemGNFS> smooths;
  long numSmooths = 0, rationalPrimes = 0, algebraicPrimes = 0;
  uint64_t entries = 0ull;
  entries = sieve(n, m, f, logMaxP2, L, b, maxA, logm, t, numSmooths, rationalPrimes, algebraicPrimes, primes, RFB, AFB, smooths)
                    .first;

  uint32_t nCols = numSmooths, nRows = rationalPrimes + algebraicPrimes + t;
  std::vector<uint32_t> mat;
  mat.reserve(entries<<1);
  uint64_t entries2 = getMatrix(mat, smooths, rationalPrimes, algebraicPrimes, QCB);
  std::vector<uint64_t> result;
  result.resize(numSmooths);
  uint32_t Nsol = blanczos(mat.data(), entries2, nRows, nCols, result.data());

  uint64_t totSol = (1ull<<(Nsol -1ull)) - 1ull;
  totSol = (totSol<<1ull) + 1ull;
  if(Nsol == 0) totSol = 0;

  #ifdef DEBUG
  if(totSol > 33u) totSol = 33u;
  #endif
  std::vector<uint32_t> U;
  for(uint64_t mask = 1ull; mask < totSol; mask++){
    std::cout << mask << std::endl;
    for(uint32_t i = 0u; i < numSmooths; i++){
        bool taken = ((std::bitset<64> (mask & result[i]).count() % 2) == 1);
        if(taken){
            U.push_back(i);
        }
    }

    // controllare che il prodotto di a-bm sia un quadrato
    factorization ratFat;
    thread_local ZZ P, num, tmp;
    bool algebraic = false;
    for(uint32_t i : U){
        if(!algebraic){
            num = smooths[i].a - smooths[i].b * m;
        } else{
            num=0;
            tmp=1;
            for(long j=0; j<=deg(f); j++){
                num += coeff(f, j)*tmp*power(smooths[i].b, deg(f)-j);
                tmp*=smooths[i].a;      // lo moltiplico con un accumulatore perché potrebbe essere 0
            }
        }
        if(num < 0){
          num = abs(num);
          ratFat[conv<ZZ>(-1)]++;
        }

        for(auto& primo : primes){
            P = primo.first;
            while(num%P == 0){
                ratFat[P]++;
                num/=P;
            }
        }

        if(num > 1) ratFat[num]++;
    }

    for(auto& coso : ratFat){
        ASSERT_TRUE((coso.second % 2) == 0) << "Prime p=" << coso.first << " has odd degree in the product of a-bm";
    }


    // controllare che il prodotto delle norme sia un quadrato
    factorization algFat;
    algebraic = true;
    for(uint32_t i : U){
        if(!algebraic){
            num = smooths[i].a - smooths[i].b * m;
        } else{
            num=0;
            tmp=1;
            for(long j=0; j<=deg(f); j++){
                num += coeff(f, j)*tmp*power(smooths[i].b, deg(f)-j);
                tmp*=smooths[i].a;      // lo moltiplico con un accumulatore perché potrebbe essere 0
            }
        }
        if(num < 0){
          num = abs(num);
          algFat[conv<ZZ>(-1)]++;
        }

        for(auto& primo : primes){
            P = primo.first;
            while(num%P == 0){
                algFat[P]++;
                num/=P;
            }
        }

        if(num > 1) algFat[num]++;
    }

    for(auto& coso : algFat){
        //std::cout<<coso.first<<" "<<coso.second<<" yee\n";
        EXPECT_TRUE((coso.second % 2) == 0) << "Prime p=" << coso.first << " has odd degree in the product of N(a, b)";
    }


    // dovrei controllare anche i quadratic characters...

  }

}


TEST(GNFS, IPB){
  ZZ n(2601699059);       // questo n genera un polinomio f con radici multiple modulo almeno 4 primi p, di cui uno dà esponente 4
  mpz_class nMPZ = 2601699059;

  long d, k, l, t;
  ZZ m, B;
  ZZX f;
  chooseParams(nMPZ, d, m, f, B);

  ZZ approxSizeX = conv<ZZ>("14272809293693672760236472311799071135480407340249456409739060630501718744859577807931963255944303817439600981719075140666405677019619267901277870858178618119743873407633109288774679593737780130609434906310650203342053608298307114931428407828641589656896990339587685062445103299971959988673570146435627636662672823036279746539829357096909005032932471318588477554189580069502760522665091841972483256163067600440067555756560726271436226408770502440639815161267283733103462972385810397545849300000933799258565095946576701402211193850929303677608722460439167668820423613379887943909376");
  ZZ prodInIPB(1);
  ZZ last(0), inertPrime;
  std::vector<ZZ> IPB;
  ZZ_pX fp;

  while(prodInIPB < approxSizeX){
    findNextInertPrime(f, last, inertPrime);

    IPB.push_back(inertPrime);
    prodInIPB*=inertPrime;
    last = inertPrime;
  }

  for(auto& P : IPB){
    // verificare che f sia irriducibile mod p
    ZZ_p::init(P);
    fp = conv<ZZ_pX>(f);

    ASSERT_TRUE(ProbIrredTest(fp, 20)) << "Prime p=" << P << " in IPB is not inert";
  }

}




#ifndef DEBUG
TEST(GNFS, bigInput) {
  mpz_class n;
  n = "2694510496740314556501037319858087";   // 100710366161784104467949462748056679263
  mpz_class p, q;
  gnfs(n, p, q);

  ASSERT_EQ(n, p*q);
}
#endif