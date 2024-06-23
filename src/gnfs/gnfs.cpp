#include "gnfs.hpp"

#ifdef DEBUG
#include "../utils/utils_print/print_stuff.hpp"
#include <cassert>
#include<iostream>
#endif
#include<iostream>
#include<sys/resource.h>

#include "../utils/utils_mpz/mpz_ZZ.hpp"
#include "../utils/utils_NTL/utils_ZZX.hpp"
#include "factorBases.hpp"



void gnfs(const mpz_class& nMPZ, mpz_class& fattore1, mpz_class& fattore2){
    // assume n>216, n dispari, n composto, n non potenza


    // 1. preparazione

    // 1.a convertire n in uno ZZ
    ZZ n;
    mpz_2_ZZ(nMPZ, n);
    std::cout << n << std::endl;

    /*mpz_class coso;
    mpz_from_ZZ(n, coso);
    std::cout << coso << std::endl;*/   // la conversione sembra funzionare

    // 1.b scelta parametri: costruire polinomio (e testare irriducibilità)
    long d;
    ZZ m, B;
    ZZX f;
    chooseParams(nMPZ, d, m, f, B);

    #ifdef DEBUG
    std:: cout << "m: " << m << std::endl;
    std::cout << "f: " << f << std::endl;
    {
        ZZ coso;
        ZZX_eval(coso, f, m);
        std::cout << "f(m)==n? " << ((coso==n) ? "true" : "false") <<std::endl;
    }
    #endif

    ZZX fPrime;
    ZZ tmp(0);
    // provare a fattorizzare f
    // controllare gcd(m, n) e gcd(f'(m), n)
    // calcolare f'
    bool foundEarlyFactor = findEarlyFactors(n, tmp, f, fPrime, m);
    #ifdef DEBUG
    if(foundEarlyFactor)
        std:: cout << "Found early factors of n" << std::endl;
    else
        std:: cout << "Did not find early factors of n" << std::endl;
    #endif
    if(foundEarlyFactor){
        mpz_from_ZZ(tmp, fattore1);
        fattore1 = abs(fattore1);
        fattore2 = nMPZ/fattore1;
        return;
    }


    // 1.c costruire factor bases
    factorBase RFB, AFB, QCB;
    ZZ L;                                                           // large prime bound, forse serve solo dentro buildFactorBases
    uint8_t logMaxP2 = buildFactorBases(n, f, RFB, AFB, QCB, B, L);    // log of large prime bound



    // 2. sieving

    // 2.a sieving razionale

    // 2.b sieving algebrico




    // 3. dipendenze lineari

    // 3.a preprocessing (trial division, segno, e quadratic characters)

    // 3.b block lanczos




    // iterare sulle dipendenze lineari trovate

    // 4. radici

    // 4.a radice razionale (includendo f'(m))

    // 4.b trovare finite fields applicabili (volendo si può fare durante costruzione factor bases)/aggiungerne altri se non bastano

    // 4.c radice in finite fields, metodo di Couveignes

    // 4.d CRT




    // 5. fattorizzare n

    // 5.a GCD

    // 5.b convertire p e q in mpz_class per metterli in fattore1 e fattore2




    // 5. gestire fail...
}








uint8_t buildFactorBases(const ZZ& n, const ZZX& f, factorBase& RFB, factorBase& AFB, factorBase& QCB, const ZZ& B, ZZ& L){
    std::vector<std::pair<long, uint8_t>> primes;
    uint8_t logMaxP = genPrimesList(primes, B);
    L = primes.back().first;
    L*=L;

    // build RFB
    // pol g = x - m    ((a/b)-m)modp == (a-bm)modp (b!=0modp)
    // primi fino a B, includere potenze fino a 2B


    // build AFB
    // pol f
    // primi fino a B, includere potenze fino a 2B



    // build QCB
    // pol f, assicurarsi di non beccare radici multiple
    // solo q>L, e bisogna generarne 3*log2(n)  (o meno?)
    long t = 3*log2(conv<double>(n));
    buildQCB(QCB, f, L, t);


    return (logMaxP << 1);      // log di L
}






bool findEarlyFactors(const ZZ& n, ZZ& fattore, const ZZX& f, ZZX fPrime, const ZZ& m){
    // provare a fattorizzare f
    // controllare gcd(m, n) e gcd(f'(m), n)
    // calcolare f'

    fattore = GCD(m, n);    // m non può essere +-n
    if(fattore > 1) return true;

    diff(fPrime, f);
    ZZ tmp;
    ZZX_eval(tmp, fPrime, m);
    fattore = GCD(tmp, n);
    if(fattore > 1) return true;

    vec_pair_ZZX_long factors;
    factor(tmp, factors, f);

    for(auto& [coso, i] : factors){     // assumo che se c'è stata una fattorizzazione non banale di f allora questa ne induce una non banale di n,
                                        // e in particolare che se non è stata trovata alcuna fattorizzazione non banale di n, allora f è irriducibile
                                        // (dovrebbe essere vero perché il termine noto di f è <= m/2)
        ZZX_eval(fattore, coso, m);
        fattore = abs(fattore);
        if(1 < fattore && fattore < n) return true;
    }

    return false;
}



void chooseParams(const mpz_class& n, long& d, ZZ& m, ZZX& f, ZZ& B){
    double eps = 1.0/128;              // epsilon > 0 arbitrario ma fisso
    double beta = cbrt(8.0/9) + eps;

    double logn = log(n.get_d()), loglogn = log(logn);

    // d = [(2/beta)^(1/2)*(logn/loglogn)^(1/3)], d dev'essere >= 3 e dispari       // grado del polinomio
    double deg = sqrt(2/beta)*cbrt(logn/loglogn);
    //std::cout << int(deg) <<std::endl;
    d = floor(deg);
    d = max(d, 3);
    if(d%2 == 0) d++;

    // B = exp(beta*(logn)^(1/3)*(loglogn)^(2/3)))          // bound per i primi in RFB e AFB
    double bound = exp(beta*cbrt(logn*loglogn*loglogn));
    if(bound < 100) bound = 100.0;
    conv(B, bound);
    //std::cout<< B << std::endl;

    // ricerca del polinomio monico f di grado d a coefficienti più piccoli per cui f(m)=n
    mpz_class mstart, mend, mbest, mcurr, bestCoeffSize = n, currCoeffSize, lastCoeff, currCoeff, tmp;
    mpz_root(mstart.get_mpz_t(), n.get_mpz_t(), d);
    mbest = mstart;
    tmp = n/2;
    mpz_root(mend.get_mpz_t(), tmp.get_mpz_t(), d);
    mend++;

    if(mend < mstart - 1'000'000) mend = mstart - 1'000'000;    // facciamo massimo 1mln di prove

    for(mcurr = mstart; mcurr >= mend; mcurr--){
        tmp = n;
        currCoeffSize = 0;
        lastCoeff = 0;
        // estrarre coeff da c_0 a c_(d-1), poi verificare c_d==1
        for(long i=0; i<=d; i++){
            currCoeff = tmp%mcurr;
            tmp/=mcurr;

            if(i!=d && lastCoeff > (mcurr/2)){
                lastCoeff -= mcurr;
                currCoeff++;
            }

            if(abs(lastCoeff) > currCoeffSize){
                currCoeffSize = abs(lastCoeff);
            }
            lastCoeff = currCoeff;
        }

        if(lastCoeff != 1){
            mpz_pow_ui(tmp.get_mpz_t(), mcurr.get_mpz_t(), d);
            std::cerr << "Errore nella scrittura in base m " << n/tmp << " " << lastCoeff << std::endl;
            exit(3);
        }

        if(currCoeffSize < bestCoeffSize){
            bestCoeffSize = currCoeffSize;
            mbest = mcurr;
        }
    }

    mpz_2_ZZ(mbest, m);
    f.SetLength(d+1);
    ZZ tmpZZ;
    mpz_2_ZZ(n, tmpZZ);
    for(long i=0; i<=d; i++){
        SetCoeff(f, i, tmpZZ%m);
        tmpZZ /= m;

        if(0<i && i<d && coeff(f, i-1) > m/2){
            SetCoeff(f, i-1, coeff(f, i-1)-m);
            SetCoeff(f, i, coeff(f, i)+1);
        }
    }

}