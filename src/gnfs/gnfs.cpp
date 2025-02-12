#include "gnfs.hpp"

#ifdef DEBUG
#include "../utils/utils_print/print_stuff.hpp"
#include <cassert>
#include<iostream>
#endif
#include<iostream>
//#include<sys/resource.h>

#include "../utils/statsVar.hpp"
#include "../utils/utils_mpz/mpz_ZZ.hpp"
#include "../utils/utils_NTL/utils_ZZX.hpp"
#include "../utils/utils_NTL/hash_ZZ.hpp"
#include "factorBases.hpp"
#include "sieving.hpp"
#include "computeSquareRoots.hpp"



void gnfs(const mpz_class& nMPZ, mpz_class& fattore1, mpz_class& fattore2){
    // assume n>216, n dispari, n composto, n non potenza

    timeval t1;
    if(printStats){
        std::clog << "{";
        initT1(t1);
    }


    // 1. preparazione

    // 1.a convertire n in uno ZZ
    ZZ n;
    mpz_2_ZZ(nMPZ, n);
    #ifdef DEBUG
    std::cout << n << std::endl;
    #endif

    /*mpz_class coso;
    mpz_from_ZZ(n, coso);
    std::cout << coso << std::endl;*/   // la conversione sembra funzionare

    // 1.b scelta parametri: costruire polinomio (e testare irriducibilità)
    long d;
    ZZ m, B;
    ZZX f;
    chooseParams(nMPZ, d, m, f, B);

    if(printStats){
        printElapsedTime("Time to choose parameters", t1);
        std::clog << "\"m\": " << m << ", "
                    << "\"B\": " << B << ", "
                    << "\"d\": " << d << ", "
                    << "\"f\": ["
                    ;
        for(long i = 0; i<=d; i++){
            std::clog << coeff(f, i);
            if(i < d) std::clog << ", ";
        }
        std::clog << "], ";
    }
    #ifdef DEBUG
    std::cout << "m: " << m << std::endl;
    std::cout << "f: " << f << std::endl;
    {
        ZZ coso;
        ZZX_eval(coso, f, m);
        std::cout << "f(m)==n? " << ((coso==n) ? "true" : "false") <<std::endl;
    }
    std::cout << "B: " << B << std::endl;
    #endif

    ZZX fPrime;
    ZZ tmp(0);
    // provare a fattorizzare f
    // controllare gcd(m, n) e gcd(f'(m), n)
    // calcolare f'
    bool foundEarlyFactor = findEarlyFactors(n, tmp, f, fPrime, m);
    #ifdef DEBUG
    if(foundEarlyFactor)
        std::cout << "Found early factors of n" << std::endl;
    else
        std::cout << "Did not find early factors of n" << std::endl;
    #endif
    if(foundEarlyFactor){

        mpz_from_ZZ(tmp, fattore1);
        fattore1 = abs(fattore1);
        fattore2 = nMPZ/fattore1;

        std::clog << "\"outcome\": \"Found early factors from parameters\"";
        std::clog << "}\n";

        return;
    }


    // 1.c costruire factor bases
    factorBase RFB, AFB, QCB;
    long k, l, t;                   // dimensioni di RFB, AFB e QCB senza potenze
    ZZ L;                                                           // large prime bound, forse serve solo dentro buildFactorBases
    primeList primes;
    uint8_t logMaxP2;
    std::tie(logMaxP2, foundEarlyFactor) = buildFactorBases(n, f, RFB, AFB, QCB, B, L, m, primes, k, l, t);    // log of large prime bound
    if(printStats){
        // stampare tempo per build factor base
        printElapsedTime("Time to build factor bases", t1);
        //std::cout << std::endl << t1.tv_usec << std::endl;

        std::clog   << "\"Large prime bound\": " << L << ", "
                    << "\"RFB size\": " << k << ", "
                    << "\"AFB size\": " << l << ", "
                    << "\"QCB size\": " << t << ", "
                    ;
    }

    if(foundEarlyFactor){
        mpz_from_ZZ(L, fattore1);
        fattore1 = abs(fattore1);
        fattore2 = nMPZ/fattore1;

        if(printStats){
            std::clog << "\"outcome\": \"Found early factors among small primes\"";
            std::clog << "}\n";
        }

        #ifdef DEBUG
        std::cout << "Found early factors of n among small primes" << std::endl;
        #endif

        return;
    }

    #ifdef DEBUG
    std::cout << "Number of primes: " << primes.size() << std::endl;
    std::cout << "RFB size: " << RFB.size() << ", k: " << k << std::endl;
    std::cout << "AFB size: " << AFB.size() << ", l: " << l << std::endl;
    std::cout << "QCB size: " << QCB.size() << ", t: " << t << std::endl;
    std::cout << "Large prime bound: " << L << std::endl;

    //std::cout << "AFB:" << std::endl;
    //printFBgnfs(AFB);
    #endif




    // 2. sieving

    long maxA;                                  // metto un limite per assicurarmi che sia un long e che non si ripetano valori nel sieving razionale
    if(B >= m/2) maxA = conv<long>((m-1)/2);         // |a|<=m/2 --> |a-bm| <= m(1/2+b) --> log|-a-bm| <= logm + log(1/2+b) < logm + logb + 1
    else maxA = conv<long>(B);                  //          --> |a-bm| >= m(b-1/2) --> log|a-bm| >= logm + log(b-1/2) >= logm + logb - 1
                                                // --> uso logm + logb + 1 - logMaxP2 come quantità da raggiungere nel sieving razionale
    uint8_t logm = log2(conv<double>(m));

    ZZ b(1);


// iterazioni esterne in caso di fail di blanczos
for(int iterazioniEsterne = 0; iterazioniEsterne<20; iterazioniEsterne++){

    std::vector<smoothElemGNFS> smooths;
    long numSmooths = 0, rationalPrimes = 0, algebraicPrimes = 0;
    uint64_t entries = 0ull;
    ZZ fatt1sieve(0);
    std::tie(entries, fatt1sieve) = sieve(n, m, f, logMaxP2, L, b, maxA, logm, t, numSmooths, rationalPrimes, algebraicPrimes, primes, RFB, AFB, smooths);
    if(printStats){
        // tempo per il sieving (le altre stat le ho printate in sieve)
        printElapsedTime("Time to activate sieve", t1);
    }
    if(fatt1sieve != 0){
        mpz_from_ZZ(fatt1sieve, fattore1);
        fattore1 = abs(fattore1);
        fattore2 = nMPZ/fattore1;

        if(printStats){
            std::clog << "\"outcome\": \"Found early factors among large primes\"";
            std::clog << "}\n";
        }

        #ifdef DEBUG
        std::cout << "Found early factors of n among large primes" << std::endl;
        #endif

        return;
    }


    // 3. dipendenze lineari

    uint32_t nCols = numSmooths, nRows = rationalPrimes + algebraicPrimes + t;
    std::vector<uint32_t> mat;
    mat.reserve(entries<<1);        // num minimo di entries che ci saranno in mat, ma poi ci andranno anche i quadratic characters

    // 3.a preprocessing (trial division (già fatta in sieve), segno (già fatto in sieve), e quadratic characters)
    uint64_t entries2 = getMatrix(mat, smooths, rationalPrimes, algebraicPrimes, QCB);  // dopo la chiamata, in smooths resteranno solo le coppie (a,b) per risparmiare memoria
    if(printStats){
        printElapsedTime("Time to build matrix", t1);
        std::clog   << "\"Matrix rows\": " << nRows << ", "
                    << "\"Matrix cols\": " << nCols << ", "
                    << "\"Matrix entries\": " << entries2 << ", "
                    << "\"Entries per row\": " << (entries2/double(nRows)) << ", "
        ;
    }

    // 3.b block lanczos
    std::vector<uint64_t> result;
    result.resize(numSmooths);
    uint32_t Nsol = blanczos(mat.data(), entries2, nRows, nCols, result.data());

    uint64_t totSol = (1ull<<(Nsol -1ull)) - 1ull;      // un po' sus ma dovrebbe andare
    totSol = (totSol<<1ull) + 1ull;
    if(Nsol == 0) totSol = 0;

    if(printStats){
        // stampare tempo per algebra lineare
        printElapsedTime("Time for linear algebra", t1);

        std::clog   << "\"Nullspace vectors\": " << Nsol << ", "
                    << "\"Span size:\": " << totSol << ", "
        ;
    }

    #ifdef DEBUG
    std::cout << "Dimensione base nullspace trovata: " << Nsol << std::endl;
    #endif
    if(Nsol == 0) std::cerr << "\"Warning\": \"0 solutions were found\"" << std::endl;



    ZZ_p y, x;
    ZZ_pX fPrimeMod;
    std::vector<uint32_t> U;     // indici delle coppie smooth prese
    factorization ratFat, algFat;
    std::vector<ZZ> IPB;        // inert prime base
    ZZ prodInIPB(1);
    ZZ approxSizeX, last(deg(f)), inertPrime;


    // iterare sulle dipendenze lineari trovate
    for(uint64_t mask = 1ull; mask < totSol; mask++){
        U.clear();
        ratFat.clear();
        algFat.clear();

        for(uint32_t i = 0u; i < numSmooths; i++){
            bool taken = ((std::bitset<64> (mask & result[i]).count() % 2) == 1);
            if(taken){
                U.push_back(i);
            }
        }
        if(printStats){
        std::clog   << "\"U size\": " << U.size() << ", "
        ;
    }
        #ifdef DEBUG
        std::cout << "Dimensione di U: " << U.size() << std::endl;
        #endif


        // 4. radici

        // 4.0 ottenere fattorizzazione del prodotto degli a-bm e fattorizzazione del prodotto delle norme
        // = mappe da ZZ a long

        computeFactorization(m, f, primes, U, smooths, ratFat);
        computeFactorization(m, f, primes, U, smooths, algFat, true);


        // 4.a radice razionale (includendo f'(m))
        ZZ_p::init(n);
        y = 1;
        computeProdMod(n, ratFat, y);
        y *= eval(conv<ZZ_pX>(fPrime), conv<ZZ_p>(m));
        /*if(y==0){
            std::cout << "f'(m) mod n: " << eval(conv<ZZ_pX>(fPrime), conv<ZZ_p>(m)) << std::endl;
        }*/


        // 4.b trovare finite fields applicabili /aggiungerne altri se non bastano
        estimateXSize(f, U, smooths, approxSizeX);

        // approxSizeX rischia di diventare enorme quindi non lo stampo in printStats per non intasare i file
        #ifdef DEBUG
        std::cout << "Stima dimensione di x: " << approxSizeX << std::endl;
        // 14272809293693672760236472311799071135480407340249456409739060630501718744859577807931963255944303817439600981719075140666405677019619267901277870858178618119743873407633109288774679593737780130609434906310650203342053608298307114931428407828641589656896990339587685062445103299971959988673570146435627636662672823036279746539829357096909005032932471318588477554189580069502760522665091841972483256163067600440067555756560726271436226408770502440639815161267283733103462972385810397545849300000933799258565095946576701402211193850929303677608722460439167668820423613379887943909376
        #endif
        while(prodInIPB < approxSizeX){
            findNextInertPrime(f, last, inertPrime);

            IPB.push_back(inertPrime);
            prodInIPB*=inertPrime;
            last = inertPrime;

            /*#ifdef DEBUG
            std::cout << "Found an inert prime: " << inertPrime << std::endl;
            #endif*/

            if(n%inertPrime == 0){
                mpz_from_ZZ(inertPrime, fattore1);
                fattore1 = abs(fattore1);
                fattore2 = nMPZ/fattore1;

                if(printStats){
                    std::clog << "\"outcome\": \"Found early factors among inert primes\"";
                    std::clog << "}\n";
                }

                #ifdef DEBUG
                std::cout << "Found early factors of n among inert primes" << std::endl;
                #endif

                return;
            }
        }
        if(printStats){
            // stampare tempo per costruire IPB
            printElapsedTime("Time to build IPB", t1);

            std::clog   << "\"IPB size\": " << IPB.size() << ", "
            ;
        }
        #ifdef DEBUG
        std::cout << "IPB size: " << IPB.size() << std::endl;
        #endif


        // 4.c radice in finite fields, metodo di Couveignes, 4.d e contemporaneamente CRT

        ZZ kl, Pl, bk;//, esp1, esp2;
        ZZ_pX betaL;
        //ZZ_p tmpP;

        ZZX esse(0), zeta(0);
        std::vector<double> kappa(deg(f)+1, 0);

        for(auto& l : IPB){
            ZZ_pPush push(l);       // serve per calcolare fl e poi lavorare in ZZ_pE

            Pl = prodInIPB/l;
            kl = InvMod(Pl%l, l);

            compSquareRootInGF(f, l, fPrime, algFat, U, smooths, betaL);

            for(long i=0; /*i<=deg(f) &&*/ i<=deg(betaL); i++){
                bk =(conv<ZZ>(coeff(betaL, i))*kl)%l;

                kappa[i] += conv<double>(bk)/conv<double>(l);
                SetCoeff(esse, i, (coeff(esse, i)+(bk*Pl)%n)%n);
            }

        }

        for(long i=0; i<=deg(f); i++){
             SetCoeff(zeta, i, conv<ZZ>(kappa[i] + 0.5));   // aggiunto 0.5 per ottenere rounding invece di floor
        }

        ZZ_p::init(n);
        ZZ_pX pol;
        pol = conv<ZZ_pX>(esse) - conv<ZZ_p>(prodInIPB)*conv<ZZ_pX>(zeta);
        x = eval(pol, conv<ZZ_p>(m));

        #ifdef DEBUG
        std::cout << "x^2-y^2 mod n: " << x*x-y*y << std::endl;
        #endif

        // 5. fattorizzare n

        // 5.a GCD
        ZZ fatt1, fatt2;
        fatt1 = abs(GCD(conv<ZZ>(x-y), n));
        fatt2 = abs(GCD(conv<ZZ>(x+y), n));
        /*std::cout << fatt1 << " " << fatt2 << std::endl;
        if(fatt1==n && fatt2 == n){
            std::cout << "x: " << x << ", y: " << y << std::endl;
        }*/

        if(fatt1 == n) fatt1 = n/fatt2;     // potrebbero essere n e q --> vogliamo p e q
        if(fatt2 == n) fatt2 = n/fatt1;     // potrebbero essere p ed n --> vogliamo p e q

        if(fatt1 != 1 && fatt2 != 1){         // adesso se sono entrambi diversi da 1 sono anche entrambi diversi da n
        // 5.b convertire p e q in mpz_class per metterli in fattore1 e fattore2
            mpz_from_ZZ(fatt1, fattore1);
            mpz_from_ZZ(fatt2, fattore2);

            if(printStats){
                // stampare tempo per ultima fase
                printElapsedTime("Time for roots and gcd", t1);

                std::clog << "\"Non-trivial congruence position\": " << mask << ", ";
                std::clog << "\"outcome\": \"Found smooth squares\"";
                std::clog << "}\n";
            }

            #ifdef DEBUG
                std::cout << "Trovati fattori in posizione: " << mask << std::endl;
            #endif
            return;
        }
    }

    // 5. gestire fail...
    if(printStats){
        std::clog << "\"outcome " << (iterazioniEsterne+1) << "\": \"Failed\", ";
    }

    // svuotare le strutture dati usate solo dal sieving in poi,    ma in realtà non serve perché sono dichiarati dentro al for, quindi si svuotano da soli a ogni iterazione
    //smooths.clear();
    //usedPrimes.clear();
    //mat.clear();
    //result.clear();
}

    if(printStats){
        std::clog << "\"outcome\": \"Failed all 20 tries\"";
        std::clog << "}\n";
    }

}





uint64_t getMatrix(std::vector<uint32_t>& mat, std::vector<smoothElemGNFS>& smooths, long numRatPrimes, long numAlgPrimes, const factorBase& QCB){
    // dopo la chiamata, in smooths resteranno solo le coppie (a,b) per risparmiare memoria
    uint64_t entries = 0ull;

    uint32_t i = 0, offset;
    for(auto& coppia : smooths){
        // inserire fattori razionali
        offset = 0;
        for(auto& pos : coppia.rationalPpos){
            entries++;
            mat.push_back(pos + offset);        // la riga è l'indice del primo
            mat.push_back(i);                   // la colonna è l'indice del numero smooth
        }
        coppia.rationalPpos.clear();


        // inserire fattori algebrici
        offset = numRatPrimes;
        for(auto& pos : coppia.algebraicPpos){
            entries++;
            mat.push_back(pos + offset);     // la riga è l'indice del primo
            mat.push_back(i);       // la colonna è l'indice del numero smooth
        }
        coppia.algebraicPpos.clear();

        // inserire quadratic characters
        offset += numAlgPrimes;
        long j = 0;
        for(auto& qc : QCB){
            thread_local ZZ coso, coso2;
            coso = Jacobi((coppia.a - coppia.b * qc.r)%qc.p, qc.p);
            #ifdef DEBUG
            coso2 = PowerMod((coppia.a - coppia.b * qc.r)%qc.p, (qc.p-1)/2, qc.p);
            if(coso2 > 1) coso2 -= qc.p;
            //if(coso != coso2) std::cout << "coso: " << coso << ", coso2: " << coso2 << std::endl;
            assert(coso == coso2);
            #endif
            if(coso < 0){       // se il simbolo di Jacobi è -1, aggiungiamo il quadratic character nella matrice
                entries++;
                mat.push_back(j + offset);      // la riga è l'indice del primo
                mat.push_back(i);               // la colonna è l'indice del numero smooth
            }

            j++;
        }
        i++;
    }

    return entries;
}


std::pair<uint64_t, ZZ> sieve(const ZZ& n, const ZZ& m, const ZZX& f, const uint8_t logMaxP2, const ZZ& L, ZZ& b, const long maxA, const uint8_t logm, const long quadChars, long& numSmooths, long& rationalPrimes, long& algebraicPrimes,  const primeList& primes, const factorBase& RFB, const factorBase& AFB, std::vector<smoothElemGNFS>& smooths){
    uint64_t entries = 0ULL;      // serve davvero?

    #ifdef DEBUG
    std::cout << "Start sieving, initial b=" << b << std::endl;
    long modForPrint = 10;
    #endif

    initSeed();                  // uno nuovo ogni volta che dobbiamo far ripartire il setaccio

    sieveArray setaccio(2*maxA + 1, (uint8_t) 0u);
    std::vector<bool> isProbSmooth(2*maxA + 1);   // sarebbbero indicizzati da -maxA a +maxA
    std::vector<ZZ> probSmooths;

    long nbit = (long) (log2(conv<double>(n))),
        threshold = std::max((long) 64, nbit>>1),       // per blanczos ne servono almeno 64 in più
        newRatPrimes = 0, newAlgPrimes = 0, partialSmooths = 0;             // solo per statistiche

    rationalPrimes = 1; algebraicPrimes = 1;        // sto contando -1 e -1
    //minRows = 1 + 1 + t,                            // il numero di righe sarà 1 (segno a-bm) + 1 (segno) + QCB.size() + num primi razionali usati + num ideali primi usati
    numSmooths = 0;

    //unsigned short lastLog = 0;           // può aiutare?
    ZZ maxNewP(0),                             // per statistiche (?)
    maxP(primes.back().first),
    aMinusBm, normAB, tmp;
    ZZX fb;
    long iterazioni = 0, newFound;                    // per statistiche
    double smoothsOverReportsAvg = 0;

    idealMap usedRatPrimes, usedAlgPrimes;         // comprende solo i primi usati
    usedRatPrimes[{conv<ZZ>(-1), conv<ZZ>(0)}] = 0u;
    usedAlgPrimes[{conv<ZZ>(-1), conv<ZZ>(0)}] = 0u;

    // iterazioni su b=1... finché non si trovano abbastanza coppie smooth
    for(; numSmooths < rationalPrimes + algebraicPrimes + quadChars + threshold; b++){
        iterazioni++;
        newFound = 0;

        tmp = 1;
        for(long i = deg(f); i>=0; i--){
            SetCoeff(fb, i, coeff(f, i)*tmp);
            tmp*=b;
        }

        #ifdef DEBUG
        if(b%modForPrint == 0){
            std::cout << "b: " << b << std::endl;
            std::cout << "numSmooths: " << numSmooths << std::endl;
            std::cout << "rationalPrimes: " << rationalPrimes << std::endl;
            std::cout << "algebraicPrimes: " << algebraicPrimes << std::endl;
            std::cout << "quadChars: " << quadChars << std::endl;
            std::cout << "threshold: " << threshold << std::endl;
        }
        #endif

        // chiamare line sieve
        long probSmoothFound = lineSieve(fb, logMaxP2, b, maxA, logm, RFB, AFB, setaccio, isProbSmooth, probSmooths);
        #ifdef DEBUG
        if(b%modForPrint == 0){
            std::cout << "b: " << b << "\t\tProb smooth found: " << probSmoothFound << std::endl;
        }
        #endif

        for(auto& a : probSmooths){
            a -= maxA;

            // scartare coppie con gcd(a, b)!=1;
            tmp = GCD(a, b);
            if(tmp != 1) continue;

            std::vector<std::pair<ZZ, ZZ>> fattoriRat, fattoriAlg;

            // trial division per a-bm
            aMinusBm = a-b*m;
            std::pair<bool, bool> isSmooth = trialDivide(a, b, aMinusBm, L, primes, fattoriRat);
            //std::cout << "hereherehere" << std::endl;
            if(!isSmooth.first){
                //smooths.pop_back();
                //fattoriRat.clear();


                /*#ifdef DEBUG
                if(b%modForPrint == 0)
                    std::cout << "b: " << b << "\t\tStuck here: not rat smooth" << std::endl;
                #endif
                //*/
                continue;
            }
            if(!isSmooth.second){
                if(n%fattoriRat.back().first == 0){
                    return std::make_pair(0, fattoriRat.back().first);
                }
            }
            bool wasPartialSmooth = !isSmooth.second;

            // trial division per N(a+b*theta)
            #ifdef DEBUG
            if(IsZero(normAB)) std::cout << "Norm of (" << a << ", " << b << ") is zero\nf is " << f << std::endl;
            #endif
            ZZX_eval(normAB, fb, a);
            isSmooth = trialDivide(a, b, normAB, L, primes, fattoriAlg);
            //std::cout << "hereherehere2" << std::endl;
            if(!isSmooth.first){
                //smooths.pop_back();
                //fattoriRat.clear();
                //fattoriAlg.clear();

                /*#ifdef DEBUG
                if(b%modForPrint == 0)
                    std::cout << "b: " << b << "\t\tStuck here: not alg smooth" << std::endl;
                #endif
                //*/
                continue;
            }
            if(!isSmooth.second){
                if(n%fattoriAlg.back().first == 0){
                    return std::make_pair(0, fattoriAlg.back().first);
                }
            }

            smooths.push_back({a, b});
            if(wasPartialSmooth || !isSmooth.second) partialSmooths++;      // non distinguo quelli parziali su RFB da quelli parziali su AFB

            // aggiungiamo effettivamente i fattori nelle hashmap
            std::vector<uint32_t>& indici1 = smooths.back().rationalPpos;
            moveFactors(indici1, fattoriRat, maxP, maxNewP, usedRatPrimes, newRatPrimes);
            entries+=indici1.size();

            /*
            #ifdef DEBUG
            if(b%(modForPrint*10) == 0)
                std::cout << "b: " << b << "\there 1" << std::endl;
            #endif
            //*/

            std::vector<uint32_t>& indici2 = smooths.back().algebraicPpos;
            moveFactors(indici2, fattoriAlg, maxP, maxNewP, usedAlgPrimes, newAlgPrimes);
            entries+=indici2.size();

            /*
            #ifdef DEBUG
            if(b%(modForPrint*10) == 0)
                std::cout << "b: " << b << "\there 2" << std::endl;
            #endif
            //*/

            numSmooths++;
            newFound++;
        }

        rationalPrimes = usedRatPrimes.size();
        algebraicPrimes = usedAlgPrimes.size();
        probSmooths.clear();

        if(probSmoothFound > 0){
            smoothsOverReportsAvg += newFound/((double) probSmoothFound);
        }
    }

    smoothsOverReportsAvg /= std::max(iterazioni, 1l);

    #ifdef DEBUG
    std::cerr << std::endl;
    std::cerr << "After " << iterazioni << " iteration(s):" << std::endl;
    //printSmooths(smooths);
    std::cout << "numSmooths: " << numSmooths << std::endl;
    std::cout << "rationalPrimes: " << rationalPrimes << std::endl;
    std::cout << "algebraicPrimes: " << algebraicPrimes << std::endl;
    std::cout << "quadChars: " << quadChars << std::endl;
    std::cerr << std::endl;
    #endif

    if(printStats){
        // stampare tempo per setaccio: fuori

        std::clog   << "\"Calls to lineSieve()\": " << iterazioni << ", "
                    << "\"Num of large rational primes\": " << newRatPrimes << ", "
                    << "\"Num of large algebraic primes\": " << newAlgPrimes << ", "
                    << "\"Total rational primes\": " << rationalPrimes << ", "
                    << "\"Total algebraic primes\": " << algebraicPrimes << ", "
                    << "\"Num of true smooths\": " << (numSmooths - partialSmooths) << ", "
                    << "\"Num of partial smooths\": " << partialSmooths << ", "
                    << "\"Biggest large prime\": " << maxNewP << ", "
                    << "\"Total dim of sieving interval\": " << (iterazioni * (2*maxA + 1)) << ", "
                    ;
        ;
    }

   return std::make_pair(entries, conv<ZZ>(0));
}



std::pair<uint8_t, bool> buildFactorBases(const ZZ& n, const ZZX& f, factorBase& RFB, factorBase& AFB, factorBase& QCB, const ZZ& B, ZZ& L, ZZ& m, primeList& primes, long& k, long& l, long& t){
    uint8_t logMaxP;
    bool foundFactor = false;
    std::tie(logMaxP, foundFactor) = genPrimesList(n, primes, B);
    L = primes.back().first;
    if(foundFactor) return std::make_pair(logMaxP, foundFactor);
    L*=L;

    // build RFB
    // pol g = x - m    ((a/b)-m)modp == (a-bm)modp (b!=0modp)
    // primi fino a B, includere potenze fino a 2B
    ZZX g; SetX(g); g-=m;
    k = buildFBase(RFB, g, B, primes);

    // build AFB
    // pol f
    // primi fino a B, includere potenze fino a 2B
    l = buildFBase(AFB, f, B, primes);

    //std::cout << "built AFB" << std::endl;


    // build QCB
    // pol f, assicurarsi di non beccare radici multiple
    // solo q>L, e bisogna generarne 3*log2(n)  (o meno?)
    long t1 = 3*log2(conv<double>(n));
    std::tie(t, foundFactor) = buildQCB(n, QCB, f, L, t1);
    if(foundFactor){
        L = QCB.back().p;
        return std::make_pair(logMaxP, foundFactor);
    }

    //std::cout << "built QCB" << std::endl;

    return std::make_pair(logMaxP << 1, false);      // log di L
}






bool findEarlyFactors(const ZZ& n, ZZ& fattore, const ZZX& f, ZZX& fPrime, const ZZ& m){
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
    if(bound < 100.0) bound = 100.0;                        // garantisce che ci siano almeno 15 primi in RFB, non dà garanzie sulle altre basi però
    conv(B, bound);
    //std::cout<< B << std::endl;

    // ricerca del polinomio monico f di grado d a coefficienti più piccoli per cui f(m)=n
    mpz_class mstart, mend, mbest, mcurr, bestCoeffSize = n, currCoeffSize, lastCoeff, currCoeff, tmp;
    mpz_root(mstart.get_mpz_t(), n.get_mpz_t(), d);
    mbest = mstart;
    tmp = n/2;
    mpz_root(mend.get_mpz_t(), tmp.get_mpz_t(), d);
    mend++;

    if(mend < mstart - 10'000*((int)logn)) mend = mstart - 10'000*((int)logn);    // limito il numero di prove

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

        #ifdef DEBUG
        if(lastCoeff != 1){
            mpz_pow_ui(tmp.get_mpz_t(), mcurr.get_mpz_t(), d);
            std::cerr << "\"Error\": \"Errore nella scrittura in base m " << n/tmp << " " << lastCoeff  << "\"" << std::endl;
            exit(3);
        }
        #endif

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

    if(coeff(f, d) != 1){
        tmpZZ = power(m, d);
        ZZ nZZ; mpz_2_ZZ(n, nZZ);
        std::clog << "\"Error\": \"Errore nella scrittura in base m " << nZZ/tmpZZ << " " << coeff(f, d)  << "\"}\n" << std::endl;
        exit(3);
    }

}