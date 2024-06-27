#include "computeSquareRoots.hpp"

#include <NTL/ZZ_pE.h>
#include <algorithm>

#ifdef DEBUG
#include <cassert>
#endif

void computeFactorization(const ZZ& m, const ZZX& f, const primeList& primes, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, factorization& fattorizzazione, bool algebraic){    // algebraic è false di default
    thread_local ZZ P, num, tmp;
    for(uint32_t i : U){
        if(!algebraic){
            num = smooths[i].a - smooths[i].b * m;
        } else{
            tmp = 1; num = 0;
            for(long j=0; j<=deg(f); j++){
                num += coeff(f, j)*tmp*power(smooths[i].b, deg(f)-j);
                tmp*=smooths[i].a;      // lo moltiplico con un accumulatore perché potrebbe essere 0
            }

            #ifdef DEBUG
            thread_local ZZX aMinusBtheta;
            aMinusBtheta = 0;
            SetCoeff(aMinusBtheta, 0, smooths[i].a);
            SetCoeff(aMinusBtheta, 1, -smooths[i].b);
            assert(num == NormMod(aMinusBtheta, f));
            #endif
        }
        num = abs(num);

        for(auto& primo : primes){
            P = primo.first;
            while(num%P == 0){
                fattorizzazione[P]++;
                num/=P;
            }
        }

        if(num > 1){
            fattorizzazione[num]++;
            #ifdef DEBUG
            assert(ProbPrime(num));
            #endif
            //std::cout << num << " ";
        }

    }

    for(auto& coso : fattorizzazione){
        #ifdef DEBUG
        if(coso.second % 2)
            std::cout << std::endl << coso.first << std::endl;
        assert(coso.second%2 == 0);
        #endif
        coso.second /= 2;
    }
}


void computeProdMod(const ZZ& modulo, const factorization& fat, ZZ_p& res){
    ZZ_pPush push(modulo);
    //ZZ_p::init(modulo);
    res = 1;

    for(auto& fattore : fat){
        res *= power(conv<ZZ_p>(fattore.first), fattore.second);
        if(fattore.first%modulo == 0){
            std::cout << "Warning: " << fattore.first << " in factorization is multiple of modulus " << modulo <<std::endl;
        }
        if(modulo % fattore.first == 0){
            std::cout << "Warning: modulo=" << modulo << " is multiple of factor=" << fattore.first <<std::endl;
        }
    }
}


void estimateXSize(const ZZX& f, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, ZZ& res){
    res = 1;
    thread_local ZZ maxC, maxAB;
    maxC = 0; maxAB = 0;

    for(long i=0; i<=deg(f); i++){
        if(maxC < abs(coeff(f, i))) maxC = abs(coeff(f, i));
    }

    for(long i : U){
        if(maxAB < abs(smooths[i].a)) maxAB = abs(smooths[i].a);
        if(maxAB < abs(smooths[i].b)) maxAB = abs(smooths[i].b);
    }

    res = power(maxC*maxAB, U.size()/2);
}


void findNextInertPrime(const ZZX& f, const ZZ& last, ZZ& inertPrime){
    long d = deg(f);
    std::vector<long> primeFactsD;
    // mini-trial division per trovare i divisori primi di d, che tanto è 3 o 5 e quindi è solo d stesso, ma giusto in caso
    for(long i = 2; i<=d; i++){
        if(d%i==0){
            bool skip = false;
            for(long j : primeFactsD){
                if(i%j==0){
                    skip=true;
                    break;
                }
            }
            if(skip) continue;
            primeFactsD.push_back(i);
        }
    }

    thread_local ZZ_pX fp, g, g2;
    thread_local ZZ_pXModulus Fp;
    thread_local ZZ esp;

    inertPrime = last;
    bool found = false;
    while(!found){
        inertPrime = NextPrime(inertPrime+1);

        ZZ_pPush push(inertPrime);
        //ZZ_p::init(inertPrime);
        fp = conv<ZZ_pX>(f);
        build(Fp, fp);

        // verificare f|x^(p^d)-x   (ovvero coso mod f = 0)
        esp = power(inertPrime, d);
        g = PowerXMod(esp, Fp);             // g = x^(p^d) mod f
        SetCoeff(g, 1, coeff(g, 1)-1);      // g = x^(p^d) -x mod f
        if(!IsZero(g % fp)) continue;

        // verificare per per ogni q primo che divide d, gcd(f, x^(p^(d/q))-x)=1
        bool skip = false;
        for(long q : primeFactsD){
            esp = power(inertPrime, d/q);
            g = PowerXMod(esp, Fp);             // g = x^(p^(d/q)) mod f
            SetCoeff(g, 1, coeff(g, 1)-1);      // g = x^(p^(d/q)) -x mod f

            if(IsZero(g)){
                skip = true;
                break;
            }

            g2 = GCD(g, fp);
            MakeMonic(g2);

            if(!IsOne(g2)){
                skip = true;
                break;
            }
        }
        if(skip) continue;

        found = true;
    }
}



void compSquareRootInGF(const ZZX& f, const ZZ& l, const ZZX& fPrime, const factorization& algFat, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, ZZ_pX& res){
    thread_local ZZ esp1, esp2;
    thread_local ZZ_p tmpNl;
    thread_local ZZ_pX fl, aMinusBtheta;
    thread_local ZZ_pE fPrimeInGF, gammaL, betaL, Nl_GF;

    esp1 = 1;
    for(long i=1; i<deg(f); i++){   // (l^d-1)/(l-1) = 1 + l + l^2 + ... + l^(d-1)
        esp1 = esp1*l + 1;
    }
    esp2 = (esp1+1)/2;
    #ifdef DEBUG
    assert((esp1 + 1)%2 == 0);
    #endif

    ZZ_pPush push(l);       // serve per calcolare fl e poi lavorare in ZZ_pE
    fl = conv<ZZ_pX>(f);
    ZZ_pEPush push2(fl);

    fPrimeInGF = conv<ZZ_pE>(conv<ZZ_pX>(fPrime));
    //std::cout << l << " " << fPrime << " " << conv<ZZ_pX>(fPrime) << std::endl;

    computeProdMod(l, algFat, tmpNl);
    Nl_GF = conv<ZZ_pE>(tmpNl) * power(fPrimeInGF, esp1);       // Nl *= f'(alhpa)^((l^d-1)/(l-1))
    //Nl_GF *= conv<ZZ_pE>(power(conv<ZZ_p>(-1), (deg(f)*(deg(f)-1))/2));
    #ifdef DEBUG
    //std::cout << tmpNl << " " << fPrimeInGF << std::endl;
    res = conv<ZZ_pX>(tmpNl);
    assert(deg(res)==0);
    res = conv<ZZ_pX>(Nl_GF);
    assert(deg(res)==0);
    #endif

    gammaL = 1;
    for(long i : U){
        aMinusBtheta = 0;
        SetCoeff(aMinusBtheta, 0, conv<ZZ_p>(smooths[i].a));
        SetCoeff(aMinusBtheta, 1, -conv<ZZ_p>(smooths[i].b));

        gammaL *= conv<ZZ_pE>(aMinusBtheta);
    }
    gammaL *= (fPrimeInGF*fPrimeInGF);

    betaL = power(gammaL, esp2);
    //std::cout<<"here"<<std::endl;
    betaL /= Nl_GF;
    //std::cout<<"here2"<<std::endl;

    #ifdef DEBUG
    assert(betaL*betaL == gammaL);
    #endif

    res = conv<ZZ_pX>(betaL);
}