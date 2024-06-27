#include "factorBases.hpp"

#include <NTL/ZZ_pXFactoring.h>
#include "../utils/utils_NTL/utils_ZZX.hpp"

std::pair<uint8_t, bool> genPrimesList(const ZZ& n, primeList& primes, const ZZ& B){
    PrimeSeq s;
    long p;
    uint8_t logMaxP = 1, log2p;
    bool foundFactor = false;

    p = s.next();
    while (p <= B) {
        log2p = (uint8_t) log2(p);
        primes.push_back({p, log2p});
        logMaxP = std::max(logMaxP, log2p);

        if(n%p == 0){
            foundFactor = true;
            break;
        }

        p = s.next();
    }
    return std::make_pair(logMaxP, foundFactor);
}



void rootsOfFmodP(const ZZX& f, const ZZ p, std::vector<ZZ>& roots, bool multiple/*=true*/){
    // multiple=false serve per escludere radici multiple: quelle con fq'(s)=0 (credo)

    thread_local ZZ_p s, tmp;
    thread_local ZZ_pX fp, g, g2, der;
    thread_local ZZ_pXModulus Fp;
    thread_local ZZ s2;

    ZZ_p::init(p);              // assumo che non sia stato ancora fatto
    fp = conv<ZZ_pX>(f);
    if(!multiple){
        der = diff(fp);
    }

    // fattorizziamo gcd(fq, x^p-x) invece di fp, così ci sono solo termini di primo grado e senza molteplicità
    build(Fp, fp);
    g = PowerXMod(p, Fp);           // g = x^p mod f
    SetCoeff(g, 1, coeff(g, 1)-1);  // g = x^p -x mod f

    // qui g potrebbe essere 0, in particolare se f divide x^q-x
    if(IsZero(g)){
        g = fp;
    }

    g2 = GCD(g, fp);                // gcd(f, x^p-x) = gcd(f, x^p-x mod f)
                                    // ora g2 dovrebbe essere square-free e avere solo fattori di primo grado
    MakeMonic(g2);                  // le radici non cambiano

    Vec<ZZ_pX> factors;
    SFCanZass(factors, g2);

    tmp = 1;                        // se multiple=true, non viene mai cambiato
    for(auto& div : factors){       // div dovrebbe essere un divisore di primo grado di f
        MakeMonic(div);
        s = -coeff(div, 0);
        if(!multiple){
            tmp = eval(der, s);
        }
        if(tmp!=0){
            s2 = conv<ZZ>(s);
            s2 = s2%p;
            if(s2 < 0) s2+=p;
            roots.push_back(s2);
        }
    }
}


std::pair<long, bool> buildQCB(const ZZ& n, factorBase& QCB, const ZZX&f, const ZZ& L, long t){
    ZZ last = L, q;
    std::vector<ZZ> roots;
    bool foundFactor = false;

    long found = 0;
    while(found < t){
        q = NextPrime(last+1);
        uint8_t log2q = (uint8_t) log2(conv<double>(q));

        //ZZ_p::init(q);
        //fq = conv<ZZ_pX>(f);

        rootsOfFmodP(f, q, roots, false);

        for(auto& s : roots){
            QCB.push_back(ideal(q, s, log2q));
            found++;
        }
        roots.clear();

        last = q;
        //std::cout << q << " " << found << std::endl;

        if(n%q == 0){
            foundFactor = true;
            break;
        }
    }
    return std::make_pair(found, foundFactor);
}



long buildFBase(factorBase& FB, const ZZX&f, const ZZ& B, primeList& primes){
    long primeIdeals = 0;

    std::vector<ZZ> roots;
    ZZ tmp, Pot;
    ZZ_p invFPrime, tmp2;
    ZZ_pX fp, der;

    for(auto& primo : primes){
        ZZ p = ZZ(primo.first);
        rootsOfFmodP(f, p, roots);
        ZZ_p::init(p);
        fp = conv<ZZ_pX>(f);
        der = diff(fp);

        //if(p==5) std::cout << fp << " " << der << std::endl;

        for(auto& r : roots){
            primeIdeals++;

            // includere potenze fino a 2B (nel caso radice singola)
            tmp2 = eval(der, conv<ZZ_p>(r));

            //if(p==5) std::cout << r << " " << tmp2 << std::endl;

            if(IsZero(tmp2)){   // radice multipla, per ogni potenza p^e ne genera p^(e-1)
                //if(p==5) std::cout << "here5" << std::endl;

                // dobbiamo titrare fuori tutti i (o almeno molti) fattori p da f(r) = r^d + c_(d-1)*r^(d-1) + ... + c_1*r + c_0
                // |f(r)| <= r^d + max(|c|)*sum(r^(d-1)+...+1) < r^d + max(|c|)*2*r^(d-1) = r^(d-1) * (r + 2max(|c|)) < r^(d-1)*3max(|c|)
                // |f(r)| < p^(d-1)*3*max(|c|) < B^(d-1)*3*max(|c|) <= 3*B^(d-1)*m <= 3*B^(d-1)*n^(1/d)

                uint8_t e = 0;
                ZZX_eval(tmp, f, r);
                if(tmp == 0){
                    e = (uint8_t) 255/primo.second;
                    tmp = 1;
                }
                while(tmp%p == 0){
                    e++; tmp/=p;
                    //if(e == (uint8_t) 255) std::cout << "here" << std::endl;
                    //std::cout << "here " << e << std::endl;
                }
                //std::cout << int(e) << std::endl;

                FB.push_back(ideal(p, r, e*primo.second));      // incremento di p perché sono tutti i valori %= r mod p,
                                                                // mentre il log è e*log2p perché f(r) è diviso da p^e

            } else{ // radice singola, per ogni potenza di p ne genera solo un'altra
                invFPrime = inv(tmp2);
                Pot = p;

                while(Pot < 2*B){
                    FB.push_back(ideal(Pot, r, primo.second));
                    ZZX_eval(tmp, f, r);
                    tmp = tmp/Pot;                          // tmp = f(r)/Pot
                    tmp2 = -conv<ZZ_p>(tmp)*invFPrime;
                    r += conv<ZZ>(tmp2)*Pot;                // radice di f mod Pot*p:   r + Pot*((-f(r)/Pot * (f'(r))^(-1))%p)
                    Pot*=p;
                }
            }


        }

        roots.clear();
    }

    return primeIdeals;
}