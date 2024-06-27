#include "sieving.hpp"
#include <cstdlib>
#include <ctime>
#include "../utils/utils_NTL/utils_ZZX.hpp"

void initSeed(){
    srand(time(NULL));
}


long lineSieve(const ZZX& f, const uint8_t logMaxP2, const ZZ& b, const long maxA, const uint8_t logm, const factorBase& RFB, const factorBase& AFB, sieveArray& setaccio, std::vector<bool>& isProbSmooth, std::vector<ZZ>& probSmoothPairs){
    // assumo che il sieveArray sia lungo 2*maxA+1 e che sia tutto a zero, e stessa cosa per isProbSmooth
    long found = 0, a, p, l = 2*maxA+1;
    thread_local ZZ tmp, rad;

    // 2.a sieving razionale

    // threshold: logm + log_2(b) + 1 - logMaxP2
    uint8_t threshold; uint16_t tmpThresh;
    tmpThresh = (uint16_t) (logm + log2(conv<double>(b)) + 1);
    if(tmpThresh > logMaxP2) threshold = (uint8_t) (tmpThresh - logMaxP2);
    else threshold = (uint8_t) 0u;

    // riempire il sieveArray
    for(auto& ideale : RFB){
        // individuare radice più vicina a -maxA
        rad = (b*ideale.r)%ideale.p;
        tmp = (-conv<ZZ>(maxA))/ideale.p;       // tmp*p <= -maxA
        tmp = rad + tmp*ideale.p;
        if(tmp < -maxA) tmp += ideale.p;
        tmp += maxA;                            // per renderlo un indice valido per il sievArray
        a = conv<long>(tmp);

        // aggiungere ideale.logP a tutte le radici <= maxA
        p = conv<long>(ideale.p);
        for(; a < l; a += p){
            setaccio[a] += ideale.logP;
        }
    }

    // cercare elementi smooth
    for(a = 0; a < l; a++){
        if(setaccio[a] >= threshold){
            isProbSmooth[a] = true;
            #ifdef DEBUG
            found++;
            #endif
        }
        setaccio[a] = 0;
    }

    #ifdef DEBUG
    if(b%10 == 0){
        std::cout << "b: " << b << "\t\tProbably rational smooth: " << found << std::endl;
    }
    found = 0;
    #endif

    // 2.b sieving algebrico

    // threshold: estraggo 10 valori di a a caso tra -maxA e +maxA
    // e faccio la media aritmetica dei log_2(|fb(a)|)-logMaxP2 (che è quindi il log_2 della media geometrica dei |fb(a)|/maxP^2)
    tmpThresh = (uint16_t) 0u;
    for(long i=0; i<10; i++){
        a = rand() % l;
        tmp = 0;
        ZZX_eval(tmp, f, conv<ZZ>(a)-maxA);
        tmpThresh += (int16_t) (log2(conv<double>(abs(tmp))));
    }
    tmpThresh/=10;
    if(tmpThresh > logMaxP2) threshold = (uint8_t) (tmpThresh - logMaxP2);
    else threshold = (uint8_t) 0u;

    // riempire il sievArray
    for(auto& ideale : AFB){
        // individuare radice più vicina a -maxA
        rad = (b*ideale.r)%ideale.p;
        tmp = (-conv<ZZ>(maxA))/ideale.p;       // tmp*p <= -maxA   && (tmp+1)*p > -maxA
        tmp = rad + tmp*ideale.p;
        if(tmp < -maxA) tmp += ideale.p;
        tmp += maxA;                            // per renderlo un indice valido per il sievArray (ovvero tra 0 e 2*maxA inclusi)
        a = conv<long>(tmp);

        // aggiungere ideale.logP a tutte le radici <= maxA
        p = conv<long>(ideale.p);
        for(; a < l; a += p){
            setaccio[a] += ideale.logP;
        }
    }

    // cercare elementi smooth
    for(a = 0; a < l; a++){
        if(isProbSmooth[a] && setaccio[a] >= threshold){
            probSmoothPairs.push_back(conv<ZZ>(a));         // sto inserendo un valore tra 0 e 2*maxA inclusi, il chiamante dovrà poi fare -maxA
            found++;
        }
        setaccio[a] = 0;
        isProbSmooth[a] = false;
    }

    return found;       // non davvero utile perché poi il chiamante riesamina le coppie
}


std::pair<bool, bool> trialDivide(const ZZ& a, const ZZ &b, const ZZ& num, const ZZ& largePrimeBound, const primeList& primes, std::vector<std::pair<ZZ, ZZ>>& fattori){
    //  dice se l'elemento era smooth (anche partialmente) o no e totalmente smooth o parziale
    fattori.clear();

    thread_local ZZ tmp, P, r;
    tmp = num;

    if(tmp < 0){       // assumo che -1 sia già tra i "primi" usati
        fattori.push_back({conv<ZZ>(-1), conv<ZZ>(0)});
        tmp = -tmp;
    }

    for(auto& [p, _] : primes){
        //p = primo.first;
        P = conv<ZZ>(p);
        if(tmp%P != 0) continue;

        /*#ifdef DEBUG
        if(b%P == 0){
            std::cout << a << " " << b << " " << tmp << std::endl;
        }
        #endif*/
        r = (a*InvMod(b%P, P))%P;        // non dovrebbe dare errore perché abbiamo escluso le coppie con a e b non coprimi
        bool e = false;
        while(tmp%P == 0){
            tmp/=P;
            e = !e;
        }
        if(e){
            fattori.push_back({P, r});
        }
    }

    if(tmp > largePrimeBound){
        return std::make_pair(false, false);
    }
    if(tmp == 1){
        return std::make_pair(true, true);
    }
    // tmp dev'essere primo perché è >maxP e <(maxP)^2 e non divisibile da alcun P fino a maxP
    r = (a*InvMod(b%tmp, tmp))%tmp;
    fattori.push_back({tmp, r});

    return std::make_pair(true, false);
}


void moveFactors(std::vector<uint32_t>& indici, std::vector<std::pair<ZZ, ZZ>>& fattori, const ZZ& maxP, ZZ& maxNewP, idealMap& usedPrimes, long& newPrimes){
    indici.clear();
    for(const auto& fattore : fattori){
        if(usedPrimes.count(fattore) > 0){
            indici.push_back(usedPrimes[fattore]);
        } else{
            if(fattore.first > maxP){
                newPrimes++;
                if(fattore.first > maxNewP) maxNewP = fattore.first;
            }
            uint32_t s = usedPrimes.size();
            usedPrimes[fattore] = s;
            indici.push_back(s);
        }
    }
}


void moveFactorsForTesting(std::vector<uint32_t>& indici, std::vector<std::pair<ZZ, ZZ>>& fattori, const ZZ& maxP, ZZ& maxNewP, idealMap& usedPrimes, long& newPrimes, std::vector<ZZ> & pos2primes){
    for(auto& fattore : fattori){
        if(usedPrimes.count(fattore) > 0){
            indici.push_back(usedPrimes[fattore]);
        } else{
            if(fattore.first > maxP){
                newPrimes++;
                if(fattore.first > maxNewP) maxNewP = fattore.first;
            }
            uint32_t s = usedPrimes.size();
            usedPrimes[fattore] = s;
            indici.push_back(s);
            pos2primes.push_back(fattore.first);
        }
    }
}


uint64_t sieveForTesting(const ZZ& n, const ZZ& m, const ZZX& f, const uint8_t logMaxP2, const ZZ& L, ZZ& b, const long maxA, const uint8_t logm, const long quadChars, long& numSmooths, long& rationalPrimes, long& algebraicPrimes,  const primeList& primes, const factorBase& RFB, const factorBase& AFB, std::vector<smoothElemGNFS>& smooths, std::vector<ZZ>& pos2ratPrimes, std::vector<ZZ>& pos2algPrimes){
    uint64_t entries = 0ULL;
    initSeed();

    sieveArray setaccio(2*maxA + 1, (uint8_t) 0u);
    std::vector<bool> isProbSmooth(2*maxA + 1);
    std::vector<ZZ> probSmooths;

    long nbit = (long) (log2(conv<double>(n))),
        threshold = std::max((long) 64, nbit>>1),
        newRatPrimes = 0, newAlgPrimes = 0, partialSmooths = 0;

    rationalPrimes = 1; algebraicPrimes = 1;
    numSmooths = 0;

    ZZ maxNewP(0),
    maxP(primes.back().first),
    aMinusBm, normAB, tmp;
    ZZX fb;

    idealMap usedRatPrimes, usedAlgPrimes;
    usedRatPrimes[{conv<ZZ>(-1), conv<ZZ>(0)}] = 0u;    pos2ratPrimes.push_back(conv<ZZ>(-1));
    usedAlgPrimes[{conv<ZZ>(-1), conv<ZZ>(0)}] = 0u;    pos2algPrimes.push_back(conv<ZZ>(-1));

    for(; numSmooths < rationalPrimes + algebraicPrimes + quadChars + threshold; b++){
        tmp = 1;
        for(long i = deg(f); i>=0; i--){
            SetCoeff(fb, i, coeff(f, i)*tmp);
            tmp*=b;
        }

        lineSieve(fb, logMaxP2, b, maxA, logm, RFB, AFB, setaccio, isProbSmooth, probSmooths);

        for(auto& a : probSmooths){
            a -= maxA;

            tmp = GCD(a, b);
            if(tmp != 1) continue;

            std::vector<std::pair<ZZ, ZZ>> fattoriRat, fattoriAlg;

            aMinusBm = a-b*m;
            std::pair<bool, bool> isSmooth = trialDivide(a, b, aMinusBm, L, primes, fattoriRat);
            if(!isSmooth.first){
                continue;
            }
            bool wasPartialSmooth = !isSmooth.second;

            ZZX_eval(normAB, fb, a);
            isSmooth = trialDivide(a, b, normAB, L, primes, fattoriAlg);
            if(!isSmooth.first){
                continue;
            }

            smooths.push_back({a, b});
            if(wasPartialSmooth || !isSmooth.second) partialSmooths++;

            std::vector<uint32_t>& indici1 = smooths.back().rationalPpos;
            moveFactorsForTesting(indici1, fattoriRat, maxP, maxNewP, usedRatPrimes, newRatPrimes, pos2ratPrimes);
            entries+=indici1.size();

            std::vector<uint32_t>& indici2 = smooths.back().algebraicPpos;
            moveFactorsForTesting(indici2, fattoriAlg, maxP, maxNewP, usedAlgPrimes, newAlgPrimes, pos2algPrimes);
            entries+=indici2.size();

            numSmooths++;
        }

        rationalPrimes = usedRatPrimes.size();
        algebraicPrimes = usedAlgPrimes.size();
        probSmooths.clear();
    }

   return entries;
}