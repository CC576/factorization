#include "quadratic_sieve.hpp"

#ifdef DEBUG
#include "../utils/utils_print/print_stuff.hpp"
#include <cassert>
#include<iostream>
#endif
//#include <cassert>
//#include<algorithm>
#include<iostream>
#include<sys/resource.h>

void quadratic_sieve(mpz_class& n, mpz_class& fattore1, mpz_class& fattore2){
    struct rusage usage;
    int ret;
    timeval t1, t2, tdiff;

    if(printStats){
        std::clog << "{";
        ret = getrusage(RUSAGE_SELF, &usage);
        t1 = usage.ru_utime;
    }

    // scelta parametri
    // B e magari approx lunghezza intervallo di sieving
    mpz_class B, L;
    choose_params(n, B, L);
    #ifdef DEBUG
    std::cerr << B << " " << L << std::endl << std::endl;
    #endif
    //std::cerr<<"here0 "<< B << std::endl;

    // costruire factor base
    std::unordered_map<mpz_class, unsigned short> factorBase;
    unsigned short logMaxP2; //= buildFactorBase(n, B, factorBase);
    unsigned long long numPrimes; //= factorBase.size();

    // righe aggiunte per sperare che blanczos non si rompa; per il gnfs cambio libreria
    logMaxP2 = buildFactorBase(n, B, factorBase);
    numPrimes = factorBase.size();
    while(numPrimes < 20)
    /*do*/{
        B*=2; L*=2;
        logMaxP2 = buildFactorBase(n, B, factorBase);
        numPrimes = factorBase.size();
        //std::cerr << "numPrimes sofar: " << numPrimes << std::endl;
    } //while(numPrimes < 20);


    if(printStats){
        // stampare tempo per build factor base
        ret = getrusage(RUSAGE_SELF, &usage);
        t2 = usage.ru_utime;
        ret = timeval_subtract(&tdiff, &t2, &t1);
        std::clog   << "\"Time to build factor base\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
        t1 = t2;
        //std::cout << std::endl << t1.tv_usec << std::endl;

        std::clog   << "\"B\": " << B << ", "
                    << "\"L\": " << L << ", "
                    << "\"Factor base dim\": " << numPrimes << ", ";
    }


    #ifdef DEBUG
    std::cerr << B << " " << L << std::endl << std::endl;
    //printFactorBase(factorBase);
    #endif
    //std::cerr<<"here1"<<std::endl;

    mpz_class base = sqrt(n) + 1;   // base per la prima volta che si usa il setaccio, = primo elemento da analizzare

    // inizializzare sieve
    std::unordered_map<mpz_class, elemSetaccio> setaccio;
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n


// da qui andrebbe tutto dentro a un loop che si ripete se il gcd trova solo divisori banali di n
    for(int iterazioniEsterne = 0; iterazioniEsterne<20; iterazioniEsterne++){

    long long totRoots = initializeSieve(n, base, L, factorBase, setaccio);

    if(printStats){
        // stampare tempo per initialize sieve
        ret = getrusage(RUSAGE_SELF, &usage);
        t2 = usage.ru_utime;
        ret = timeval_subtract(&tdiff, &t2, &t1);
        std::clog   << "\"Time to initialize sieve\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
        t1 = t2;
        //std::cout << std::endl << t1.tv_usec << std::endl;

        if(iterazioniEsterne == 0){
            std::clog << "\"Max num of elements in sieve array (theoretical)\": " << totRoots << ", " ;
        }
    }

    #ifdef DEBUG
    //printSetaccio(setaccio, base);
    #endif
    //std::cerr<<"here2"<<std::endl;


    // controllare che il setaccio sia stato inizializzato correttamente
    // fatto con un test


    // attivare sieve

        // durante il setaccio i primi/le potenze vengono aggiunti a un numero dal più grande al più piccolo (mentre si analizzano i numeri precedenti)
        // --> se si inseriscono con push_front poi si ritrovano dal più piccolo al più grande (mentre si analizza quel numero)
        // --> si trovano prima i primi, poi le loro potenze; questa è un'assunzione fondamentale!!!

    std::vector<smoothElem> smooths;

    mpz_class a = -1;   // dato che nella funzione l'incremento avviene prima del resto, in questo modo il primo valore di a sarà 0
    mpz_class baseSquaredMinusN = base*base - n, twoBase = base<<1;

    unsigned int  nbit = mpz_sizeinbase(n.get_mpz_t(), 2),      // int e non long long per blanczos
                        threshold = std::max(64u, nbit>>1),     // per blanczos ne servono almeno 64 in più
                        toSearch = std::max(64u, nbit<<2),
                        numNewPrimes = 0,
                        numTotPrimes = numPrimes,               // conta tutti i primi della factor base (usati e non usati), più i large primes usati
                        missing = 0 + numTotPrimes + nbit,
                        numSmooths = 0, partialSmooths = 0;
    unsigned long long entries = 0ULL;      // serve davvero?
    unsigned short lastLog = 0;
    //#ifdef DEBUG
    unsigned short iteration = 0;
    mpz_class maxP = 0;
    //#endif

    std::unordered_map<mpz_class, uint32_t> usedPrimes;         // comprende solo i primi usati

    do{
        #ifdef DEBUG
        std::cerr << "Iteration: " << iteration << std::endl;
        std::cerr << "Threshold: " << threshold << std::endl;
        std::cerr << "toSearch: " << toSearch << std::endl;
        std::cerr << "numTotPrimes: " << numTotPrimes << std::endl;
        std::cerr << "numSmooths: " << numSmooths << std::endl;
        std::cerr << std::endl;
        //std::cout << std::endl << maxP - B*B << std::endl;
        //std::cout << std::endl << a << std::endl;
        //printSmooths(smooths);
        #endif
        //std::cerr<<"hereI"<<std::endl;

        missing = numTotPrimes + toSearch - numSmooths,         // potrebbe dare problemi questa sottrazione perché è unsigned, ma in teoria no perché toSearch >= treshold
        numSmooths += activateSieve(missing, logMaxP2, a, base, baseSquaredMinusN, twoBase, /*factorBase,*/ setaccio, smooths, lastLog);

        // preprocessing: calcolare exponent vectors, fare (o no) filtering, se non ci sono abbastanza smooth numbers tornare al sieving
            // dato ogni numero in smooths:
                // tirare via le potenze dei primi noti segnandosi quali sono quelli con esp dispari (in un qualche ordine)
                // segnarsi ciò che rimane (se diverso da 1) come nuovo primo che compare con esp 1 (quindi segnarlo nell'exp vector ridotto) e aggiungerlo ad un set di nuovi primi
            // se il numero di numeri trovati non supera numPrimes + numNewPrimes riprendere sieving, altrimenti si può ripulire il setaccio, ma segnandosi l'ultimo elem visto come nuova base
            // fare un unico set con tutti i primi che hanno degli exp dispari, e per ognuno di questi ricavare la posizione ordinata,
                // così da poter fare poi gli exponent vectors con le posizioni dei primi (notazione compatta per matrice densa)
            // non rimuovo i large primes che compaiono un'unica volta perché tanto se li tengo aggiungono sia un'incognita che un'equazione, quindi non cambia il numero di equazioni mancanti

        //elemSetaccio divisors;
        //std::cerr << numSmooths - missing << " " << numSmooths << std::endl;
        for(unsigned i = numSmooths-missing; i < numSmooths; i++){  // bisogna solo scorrere i nuovi elementi aggiunti
            unsigned int count = 0u;        // = num di fattori primi di smooths[i] con exp dispari
            auto& elem = smooths[i];
            mpz_class y = elem.y;

            auto it1 = elem.primes.before_begin(), it2 = it1;
            it2++;
            while(it2 != elem.primes.end()){
                mpz_class p = it2->first;

                unsigned int esp = mpz_remove(y.get_mpz_t(), y.get_mpz_t(), p.get_mpz_t());
                if(esp%2){
                    //divisors.push_front({p, 0});
                    //divisors.push_back({p, (short unsigned) 0});
                    if(usedPrimes.count(p) == 0){
                        uint32_t pos = usedPrimes.size();
                        usedPrimes[p] = pos;
                    }
                    entries++;
                    count++;

                    it1++; it2++;
                } else{
                    it2 = elem.primes.erase_after(it1);     // it1 resta un iteratore valido
                }
            }

            /*while(!elem.primes.empty()){
                mpz_class p = elem.primes.front().first;
                elem.primes.pop_front();
                //mpz_class p = elem.primes.back().first;
                //elem.primes.pop_back();
                unsigned int esp = mpz_remove(y.get_mpz_t(), y.get_mpz_t(), p.get_mpz_t());
                if(esp%2){
                    divisors.push_front({p, 0});
                    //divisors.push_back({p, (short unsigned) 0});
                    if(usedPrimes.count(p) == 0){
                        uint32_t pos = usedPrimes.size();
                        usedPrimes[p] = pos;
                    }
                    entries++;
                    count++;
                }
            }*/

            if(y != 1){     // y contiene esattamente un large prime
                elem.primes.push_front({y, 0});
                //divisors.push_front({y, 0});
                //divisors.push_back({y, (short unsigned) 0});
                if(usedPrimes.count(y) == 0){    // abbiamo trovato un nuovo large prime (che in un certo senso estende la factor base)
                    numNewPrimes++; // questo serve solo per statistiche in realtà
                    //std::cout << y << std::endl;
                    numTotPrimes++;
                    uint32_t pos = usedPrimes.size();
                    usedPrimes[y] = pos;
                    //#ifdef DEBUG
                    if(y>maxP) maxP=y;
                    //#endif
                }
                entries++;
                count++;
                partialSmooths++;
            }
            if(count == 0){ // y è un quadrato
                #ifdef DEBUG
                std::cerr << "Found a perfect square" << std::endl;
                #endif
                mpz_class Y = sqrt(elem.y), X = elem.x;
                fattore1 = abs(gcd(n, X-Y));                // questo pezzo deve diventare una funzione
                fattore2 = abs(gcd(n, X+Y));
                if(fattore1 == n) fattore1 = n/fattore2;   // potrebbero essere n e q --> vogliamo p e q
                if(fattore2 == n) fattore2 = n/fattore1;   // potrebbero essere p ed n --> vogliamo p e q
                if(fattore1 != 1 && fattore2 != 1){         // se sono entrambi diversi da 1 sono anche entrambi diversi da n
                    if(printStats){
                        // stampare tempo per setaccio fin qui
                        ret = getrusage(RUSAGE_SELF, &usage);
                        t2 = usage.ru_utime;
                        ret = timeval_subtract(&tdiff, &t2, &t1);
                        std::clog   << "\"Time to activate sieve\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
                        t1 = t2;
                        //std::cout << std::endl << t1.tv_usec << std::endl;

                        std::clog << "\"Dim of sieving interval\": " << (a+base - sqrt(n)) << ", ";
                        std::clog << "\"outcome\": \"Found a perfect square\", \"perfect square\": " << elem.y;
                        std::clog << "}\n";
                    }
                    return;
                }
            }
            //divisors.swap(elem.primes);     // ora elem.primes contiene solo i divisori primi con esponenete dispari,
        }                                   // che risultano in ordine crescente a parte per il large prime che viene messo all'inizio

        //#ifdef DEBUG
        iteration++;
        //#endif
    } while(numSmooths < threshold + numTotPrimes);
    setaccio.clear();   // ripulisco il setaccio per avere più spazio per l'algebra lineare

    if(printStats){
        // stampare tempo per setaccio
        ret = getrusage(RUSAGE_SELF, &usage);
        t2 = usage.ru_utime;
        ret = timeval_subtract(&tdiff, &t2, &t1);
        std::clog   << "\"Time to activate sieve\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
        t1 = t2;
        //std::cout << std::endl << t1.tv_usec << std::endl;

        std::clog   << "\"Calls to activateSieve()\": " << iteration << ", "
                    << "\"Num of large primes\": " << numNewPrimes << ", "
                    << "\"Num of true smooths\": " << (numSmooths - partialSmooths) << ", "
                    << "\"Num of partial smooths\": " << partialSmooths << ", "
                    << "\"Biggest large prime\": " << maxP << ", "
        ;
    }

    #ifdef DEBUG
    std::cerr << "After " << iteration << " iteration(s):" << std::endl;
    //printSmooths(smooths);
    std::cerr << "Threshold: " << threshold << std::endl;
        std::cerr << "toSearch: " << toSearch << std::endl;
    std::cerr << "numTotPrimes: " << numTotPrimes << std::endl;
    std::cerr << "numSmooths: " << numSmooths << std::endl;
    std::cerr << "Biggest large prime used: " << maxP << std::endl;
    std::cerr << "Number of iterations: " << iteration << std::endl;
    std::cerr << std::endl;
    //std::cout << std::endl << a << std::endl;
    #endif

    #ifdef DEBUG
    std::cerr << "Total number of primes: " << usedPrimes.size() << std::endl << std::endl;
    //printUsedPrimes(usedPrimes);
    #endif
    //std::cerr<<"here3"<<std::endl;




    // algebra lineare, block lanczos

    // la libreria restituisce fino a 64 vettori del nullspace, salvandoli in ncol interi a 64 bit
        // credo che ogni intero rappresenti i 64 valori da dare alla corrispondente incognita nelle 64 soluzioni
        // --> per iterare sulle 2^nsol soluzioni basta avere un contatore unsigned ll che va da 1 (!!!non da 0!!!) a 2^nsol-1 e fare & con ciascun intero (sepratamente)
        //      e guardare la parità del numero di bit accesi

    // nrow = usedPrimes.size(), ncol = numSmooths,
        // ovvero sulle righe ci sono i primi, sulle colonne i numeri smooth;
        // un'equazione è la somma degli esponenti mod 2 di un primo sui vari smooth,
        // una soluzione indica quali smooth usare per rendere ogni equazione 0

    std::vector<uint32_t> mat;
    mat.reserve(entries<<1);
    uint64_t entries2 = getMatrix(smooths, mat, usedPrimes);
    if(printStats){
        std::clog   << "\"Matrix rows\": " << usedPrimes.size() << ", "
                    << "\"Matrix cols\": " << numSmooths << ", "
                    << "\"Matrix entries\": " << entries2 << ", "
                    << "\"Entries per row\": " << (entries2/double(usedPrimes.size())) << ", "
        ;
    }
    #ifdef DEBUG
    std::cerr << "Number of entries: " << entries << " =? " << entries2 << std::endl;
    assert(entries == entries2);
    std::cerr << std::endl;
    assert(entries == entries2);
    #endif

    std::vector<uint64_t> result;
    result.resize(numSmooths);
    uint32_t Nsol = blanczos(mat.data(), entries, usedPrimes.size(), numSmooths, result.data());
    #ifdef DEBUG
    std::cerr << Nsol << std::endl;
    #endif
    //std::cerr<<"here4"<<std::endl;

    // differenza quadrati: radici e gcd

    uint64_t totSol = (1ull<<(Nsol -1ull)) - 1ull;      // un po' sus ma dovrebbe andare
    totSol = (totSol<<1ull) + 1ull;
    if(Nsol == 0) totSol = 0;

    if(printStats){
        // stampare tempo per algebra lineare
        ret = getrusage(RUSAGE_SELF, &usage);
        t2 = usage.ru_utime;
        ret = timeval_subtract(&tdiff, &t2, &t1);
        std::clog   << "\"Time for linear algebra\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
        t1 = t2;
        //std::cout << std::endl << t1.tv_usec << std::endl;

        std::clog   << "\"Nullspace vectors\": " << Nsol << ", "
                    << "\"Span size:\": " << totSol << ", "
        ;
    }

    #ifdef DEBUG
    std::cerr << totSol << std::endl;
    std::cerr << std::endl;
    if(Nsol == 0) std::cerr << "Warning: 0 solutions were found" << std::endl;
    #endif
    if(Nsol == 0) std::cerr << "\"Warning\": \"0 solutions were found\"" << std::endl;
    //std::cerr<<"here5" << " " << Nsol <<std::endl;

    mpz_class X, Y, Ysquared;
    for(uint64_t mask = 1ull; mask < totSol; mask++){
        X = Ysquared = 1;
        for(uint32_t i = 0u; i < numSmooths; i++){
            //uint32_t val = mask & result[i];
            //std::bitset<64> val (mask & result[i]);
            bool taken = ((std::bitset<64> (mask & result[i]).count() % 2) == 1);
            if(taken){
                X *= smooths[i].x;
                X = X%n;
                Ysquared *= smooths[i].y;
            }
        }

        Y = sqrt(Ysquared);
        fattore1 = abs(gcd(n, X-Y));                // questo pezzo deve diventare una funzione
        fattore2 = abs(gcd(n, X+Y));
        if(fattore1 == n) fattore1 = n/fattore2;   // potrebbero essere n e q --> vogliamo p e q
        if(fattore2 == n) fattore2 = n/fattore1;   // potrebbero essere p ed n --> vogliamo p e q
        if(fattore1 != 1 && fattore2 != 1){         // adesso se sono entrambi diversi da 1 sono anche entrambi diversi da n
            if(printStats){
                // stampare tempo per ultima fase
                ret = getrusage(RUSAGE_SELF, &usage);
                t2 = usage.ru_utime;
                ret = timeval_subtract(&tdiff, &t2, &t1);
                std::clog   << "\"Time for roots and gcd\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
                t1 = t2;
                //std::cout << std::endl << t1.tv_usec << std::endl;

                std::clog << "\"Dim of sieving interval\": " << (a+base - sqrt(n)) << ", ";
                std::clog << "\"Non-trivial congruence position\": " << mask << ", ";
                std::clog << "\"outcome\": \"Found product perfect square\"";
                std::clog << "}\n";
            }
            return;
        }
    }


    // in caso di fail...
    //std::cerr<<"here6"<<std::endl;
    if(printStats){
        std::clog << "\"outcome " << (iterazioniEsterne+1) << "\": \"Failed\", ";
    }

    base += a+1;                     // primo elemento da visitare in caso di fail = ultimo visitato +1
    smooths.clear();
    usedPrimes.clear();
    mat.clear();
    result.clear();
    }

    //std::cerr<<"here7"<<std::endl;
    //std::cerr<<"Something went wrong"<<std::endl;
    if(printStats){
        std::clog << "\"Dim of sieving interval\": " << (base - sqrt(n) - 1) << ", ";
        std::clog << "\"outcome\": \"Failed all 20 tries\"";
        std::clog << "}\n";
    }
}


unsigned long long getMatrix(std::vector<smoothElem>& smooths, std::vector<uint32_t>& B, std::unordered_map<mpz_class, uint32_t>& primes){   // restituisce il numero di entries
    // dopo questa funzione le liste dei fattori di ogni smooth saranno vuote
    // B deve avere già abbastanza spazio --> usare reserve

    unsigned long long entries = 0ULL;
    uint32_t n = smooths.size();


    for(uint32_t i = 0u; i < n; i++){
        auto& elem = smooths[i];
        while(!elem.primes.empty()){
            B.push_back(primes[elem.primes.front().first]);    // la riga è l'indice del primo
            //B.push_back(primes[elem.primes.back().first]);
            B.push_back(i);                        // la colonna è l'indice del numero smooth
            elem.primes.pop_front();
            //elem.primes.pop_back();

            entries++;
        }

    }
    return entries;
}



unsigned long long activateSieve(unsigned long long toFind, unsigned short maxLogp2, mpz_class& a, mpz_class& base, mpz_class& baseSquaredMinusN, mpz_class& twoBase,  /*std::unordered_map<mpz_class, unsigned short>& factorBase,*/ std::unordered_map<mpz_class, elemSetaccio>& setaccio, std::vector<smoothElem>& smooths, unsigned short& lastLog){   // il valore di ritorno è il numero di smooth trovati
    // alla fine a conterrà l'ultimo elemento ispezionato --> i setacci dopo possono ripartire da a
    thread_local mpz_class /*x,*/ y;//, P;
    unsigned long long found = 0;

    elemSetaccio divisors;

    while(found < toFind){
        a++;                    // a era l'ultimo visto
        //std::cerr<<a<<std::endl;

        auto it = setaccio.find(a);
        if(it == setaccio.end()) continue;

        // ora it->second dovrebbe essere la lista con potenze dei primi che dividono y(x(a)) e relativi logaritmi
        divisors.swap(it->second);
        setaccio.erase(it);

        // ora divisors dovrebbe essere la lista con potenze dei primi che dividono y(x(a)) e relativi logaritmi

        unsigned short sum_of_logs = 0;//, l = 0;
        for(auto& [P, l] : divisors){
            //P = coppia.first;
            //l = coppia.second;

            // passare avanti la potenza che divide y
            /*auto j = setaccio.find(a + P);
            if(j == setaccio.end()){
                //setaccio[a+P].reserve(100);
                //setaccio[a+P] = elemSetaccio{std::make_pair(P, l)};
            } else{
                j->second.push_front(std::make_pair(P, l));
            }*/
            setaccio[a+P].push_front(std::make_pair(P, l));   // se lista
            //setaccio[a+P].push_back(std::make_pair(P, l));      // se vector --> le potenze più grandi sono davanti

            // sommare i log
            sum_of_logs += l;
        }

        /*#ifdef DEBUG
        if(a%1000 == 0)
            std::cerr<< " " << (sum_of_logs>>6) << std::endl;
        #endif*/

        if(sum_of_logs < lastLog){
            divisors.clear();
            continue;
        }

        y = baseSquaredMinusN + a*(twoBase + a);
        double logyD = log2(y.get_d());
        unsigned short logy = (unsigned short) (logyD * (1<<6));

        // lastLog is < logy to allow inclusion of "large" prime    - da Wikiversity
        if(logy > maxLogp2){
            lastLog = logy - maxLogp2;
            /*#ifdef DEBUG
            std::cerr << "here " << (lastLog>>6) << " " << (maxLogp2>>6) <<std::endl;
            #endif*/
        }
        /*#ifdef DEBUG
        std::cerr<< " " << (lastLog>>6) << std::endl;
        #endif*/
        if(sum_of_logs < lastLog){
            divisors.clear();
            continue;
        }

        // abbiamo un elemento probabilmente smooth, tranne per al più un large prime
        found++;
        /*#ifdef DEBUG
        std::cerr << "found one " << found <<std::endl;
        #endif*/
        smooths.push_back(smoothElem(a+base, y));   // x = a+base, y = x^2-n
        smoothElem & elem = smooths.back();


        // durante il setaccio i primi/le potenze vengono aggiunti a un numero dal più grande al più piccolo (mentre si analizzano i numeri precedenti)
        // --> se si inseriscono con push_front poi si ritrovano dal più piccolo al più grande (mentre si analizza quel numero) (con push_back avviene al contrario)
        // --> si trovano prima i primi, poi le loro potenze; questa è un'assunzione fondamentale!!!
        // aggiustato initializeSieve per assicurarmi che ogni primo compaia prima delle potenze nelle varie liste
        /*while(!divisors.empty()){
            P = (divisors.front()).first;
            divisors.pop_front();

            if(factorBase.count(P)){    // userei contains ma c'è solo da c++20
                elem.primes.push_front(P);
            }
        }*/ divisors.swap(elem.primes);
    }

    return found;
}




long long initializeSieve(const mpz_class& n, mpz_class& base, mpz_class& L, std::unordered_map<mpz_class, unsigned short>& factorBase, std::unordered_map<mpz_class, elemSetaccio>& setaccio){
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n

    long long totRoots = 0ll;

    mpz_class p, Pot, Pot_, a, tmp1, tmp2, tmp3;
    std::vector<mpz_class> roots, roots1;
    #ifdef DEBUG
    unsigned short totLog = (short unsigned) 0;
    #endif
    for(auto& coppia : factorBase){
        p = coppia.first;
        unsigned short l = coppia.second;
        //std::cerr << p << " " << l << std::endl;
        if(p == 2){
            roots.resize(0);
            roots.push_back(1); // assumo che n sia dispari
        } else{
            roots.resize(0);
            roots.resize(2);
            a = n%p;
            findRoot(a, p, roots[0], tmp1, tmp2, tmp3);
            roots[1] = p-roots[0];
        }
        #ifdef DEBUG
        totLog += l;
        #endif

        roots1 = roots;
        //insertRoots(roots, base, p, setaccio, a, l);  // inserisco le radici mod p dopo quelle mod p^e per averle in cima alle forward list


        //potenze (ricordarsi che per p=2 c'è un'unica radice, per Pot = 4 ce ne sono 2 per quadrato, da Pot=8 in poi ogni quadrato ha  4 radici)
        Pot_ = p; Pot = Pot_ * p;
        while(Pot <= L){    // Pot è la potenza corrente, Pot_ la precedente
            a = n;
            if(p != 2){// p dispari
                // ci sono 2 radici

                extendRootP(roots[0], a, Pot_, p, tmp1, tmp2, tmp3);
                roots[1] = Pot - roots[0];

            }else if(Pot == 4){
                if(n%4 != 1) break;
                // ci sono 2 radici
                roots.resize(2);
                roots = {1, 3};
                //std::cerr<<"here"<<std::endl;
            } else if(Pot == 8){
                if(n%8 != 1) break;
                // ci sono 4 radici
                roots.resize(4);
                roots = {1, 3, 5, 7};


            } else{ // potenze di 2 successive a 8
                // ci sono 4 radici, derivano da quelle precedenti; sono r, -r, Pot_+r, Pot_-r
                extendRoot2(roots[0], a, Pot_, tmp1);
                roots[1] = Pot_ - roots[0];
                roots[2] = Pot_ + roots[0];
                roots[3] = Pot - roots[0];
            }

            totRoots += insertRoots(roots, base, Pot, setaccio, a, l);

            Pot_ = Pot;
            Pot *= p;

            #ifdef DEBUG
            totLog += l;
            #endif
        }

        totRoots += insertRoots(roots1, base, p, setaccio, a, l);
    }
    #ifdef DEBUG
    std::cerr<<"totLog: "<<totLog << " " << (totLog>>6) << std::endl <<std::endl;
    #endif
    return totRoots;
}



long long insertRoots(std::vector<mpz_class>& roots, mpz_class& base, mpz_class& P, std::unordered_map<mpz_class, elemSetaccio>& setaccio, mpz_class& a, unsigned short l){
    for(auto& r : roots){
        // individuale il valore positivo di a più piccolo per cui a+base == r mod p
        // --> a == r-base mod p
        a = (r - base) % P; // la divisione tronca verso 0
        if(a < 0) a += P;

        auto it = setaccio.find(a);
        if(it == setaccio.end()){
            setaccio[a] = elemSetaccio{std::make_pair(P, l)};
        }else{
            it->second.push_front(std::make_pair(P, l));
            //it->second.push_back(std::make_pair(P, l));
        }
    }
    return roots.size();
}



unsigned short buildFactorBase(mpz_class& n, mpz_class& B, std::unordered_map<mpz_class, unsigned short>& factorBase){
    // effettua anche controllo simbolo di Legendre, che per versione MPQS va spostato nell'inizializzazione del sieve
    mpz_class tmp;
    unsigned long long b = mpz_2_ull(B, tmp);  // sto assumendo che B sia < 2^64 e che possa stare in un long long

    std::vector<bool> composite(b+1);
    unsigned short logMaxP = 1<<6;

    for(unsigned long long p=2; p<=b; p++){
        if(composite[p]) continue;

        // p è primo
        for(unsigned long mul=p*p; mul<=b; mul+=p){
            composite[mul] = true;
        }

        mpz_from_ull(p, tmp);   // tmp contiene p

        if(p == 2 && n%2==0) continue;  // questo caso non dovrebbe verificarsi
        if(p != 2 && mpz_legendre(n.get_mpz_t(), tmp.get_mpz_t()) != 1) continue;

        double log2pD = log2(double(p));
        unsigned short log2p = (unsigned short) (log2pD * (1<<6));

        factorBase[tmp] = log2p;

        logMaxP = std::max(logMaxP, log2p);
    }
    return (logMaxP << 1); //- 0*(1<<6);     // considero large primes solo fino a (maxP^2)/2 per perdere meno tempo su numeri con fattori che non ricompariranno
}                                           // ma in realtà non dà miglioramenti significativi - cosa fare? Rimetto B^2



void choose_params(mpz_class &n, mpz_class &B, mpz_class &L){
    // B = exp((1/2 + o(1))*(logn loglogn)^(1/2)), e controllare che sia almeno 100
        // (altrimenti ci vuole troppo tempo per i numeri piccoli perché non si trovano abbastanza B-smooth)

    // L = exp((logn loglogn)^(1/2)-(1/2)*loglogn), e controllare che sia almeno 5000

    double logn = log(n.get_d()), loglogn = log(logn);
    double  b = exp((1/2.0 + 100/n.get_d())*sqrt(logn * loglogn)),                     // scelto 100/n come o(1)
            l = exp((1.0 + 10/n.get_d())*sqrt(logn * loglogn) - (1/2.0)*loglogn);     // qui usato un decimo per non prendere troppe potenze, che tanto vengono già coperte dal margine per il large prime

    B = b; L = l;

    //B*=3; L*=3;         // euristica che dà risultati molto migliori su n a 128 bit, ma rompre block lanczos su n a 160 bit...

    if(B < 100) B=100;
    if(L < 5000) L=5000;
}
