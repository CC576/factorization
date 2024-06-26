#include "computeSquareRoots.hpp"

#include <algorithm>

void computeFactorization(const ZZ& m, const ZZX& f, const primeList& primes, const std::vector<uint32_t>& U, const std::vector<smoothElemGNFS>& smooths, factorization& fattorizzazione, bool algebraic){    // algebraic è false di default
    thread_local ZZ P, num, tmp;
    for(uint32_t i : U){
        if(!algebraic){
            num = smooths[i].a - smooths[i].b * m;
        } else{
            for(long i=0; i<=deg(f); i++){
                num += coeff(f, i)*tmp*power(smooths[i].b, deg(f)-i);
                tmp*=smooths[i].a;      // lo moltiplico con un accumulatore perché potrebbe essere 0
            }
        }
        num = abs(num);

        for(auto& primo : primes){
            P = primo.first;
            while(num%P == 0){
                fattorizzazione[P]++;
                num/=P;
            }
        }

        if(num > 1) fattorizzazione[num]++;
    }

    for(auto& coso : fattorizzazione){
        coso.second /= 2;
    }
}


void computeProdMod(const ZZ& modulo, const factorization& fat, ZZ_p& res){
    ZZ_p::init(modulo);
    res = 1;

    for(auto& fattore : fat){
        res *= power(conv<ZZ_p>(fattore.first), fattore.second);
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

        ZZ_p::init(inertPrime);
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