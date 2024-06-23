#include "factorBases.hpp"

#include <NTL/ZZ_pXFactoring.h>

uint8_t genPrimesList(std::vector<std::pair<long, uint8_t>>& primes, const ZZ& B){
    PrimeSeq s;
    long p;
    uint8_t logMaxP = 1, log2p;

    p = s.next();
    while (p <= B) {
        log2p = (uint8_t) log2(p);
        primes.push_back({p, log2p});
        logMaxP = std::max(logMaxP, log2p);
        p = s.next();
    }
    return logMaxP;
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


void buildQCB(factorBase& QCB, const ZZX&f, const ZZ& L, long t){
    ZZ last = L, q;
    std::vector<ZZ> roots;

    long found = 0;
    while(found < t){
        q = NextPrime(last);
        uint8_t log2q = (uint8_t) log2(conv<double>(q));

        //ZZ_p::init(q);
        //fq = conv<ZZ_pX>(f);

        rootsOfFmodP(f, q, roots, false);

        for(auto& s : roots){
            QCB.push_back(ideal(q, s, log2q));
            found++;
        }
        roots.clear();
        roots.resize(0);

        last = q;
    }
}