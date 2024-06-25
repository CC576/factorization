#include "sieving.hpp"



long lineSieve(const ZZX& f, const uint8_t logMaxP2, const ZZ& b, const long maxA, const uint8_t logm, const factorBase& RFB, const factorBase& AFB, sieveArray& setaccio, std::vector<bool>& isProbSmooths, std::vector<ZZ>& probSmoothPairs){
    long found = 0;
    // 2.a sieving razionale
    // threshold: logm + log_2(b) + 1 - logMaxP2

    // 2.b sieving algebrico
    // threshold: estraggo 10 valori di a a caso tra -maxA e +maxA
    // e faccio la media aritmetica dei log_2(|fb(a)|)-logMaxP2 (che è quindi il log_2 della media geometrica dei |fb(a)|/maxP^2)


    return found;       // non davvero utile perché poi il chiamante riesamina le coppie
}


std::pair<bool, bool> trialDivide(const ZZ& a, const ZZ &b, const ZZ& num, const ZZ& largePrimeBound, const primeList& primes, std::vector<std::pair<ZZ, ZZ>>& fattori){
    //  dice se l'elemento era smooth (anche partialmente) o no e totalmente smooth o parziale
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

        r = (-a*InvMod(b, P))%P;        // non dovrebbe dare errore perché abbiamo escluso le coppie con a e b non coprimi
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
    r = (-a*InvMod(b, tmp))%tmp;
    fattori.push_back({tmp, r});
}


void moveFactors(std::vector<uint32_t>& indici, std::vector<std::pair<ZZ, ZZ>>& fattori, const ZZ& maxP, ZZ& maxNewP, idealMap& usedPrimes, long& newPrimes){
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
        }
    }
}