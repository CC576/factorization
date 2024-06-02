#include "quadratic_sieve.hpp"

#ifdef DEBUG
#include "../utils/utils_print/print_stuff.hpp"
#endif

void quadratic_sieve(mpz_class& n, mpz_class& fattore1, mpz_class& fattore2){

    // scelta parametri
    // B e magari approx lunghezza intervallo di sieving
    mpz_class B, L;
    choose_params(n, B, L);
    #ifdef DEBUG
    std::cerr << B << " " << L << std::endl << std::endl;
    #endif


    // costruire factor base
    std::unordered_map<mpz_class, unsigned short> factorBase;
    unsigned short logMaxP2 = buildFactorBase(n, B, factorBase);
    unsigned long long numPrimes = factorBase.size();

    #ifdef DEBUG
    printFactorBase(factorBase);
    #endif


    // inizializzare sieve
    std::unordered_map<mpz_class, elemSetaccio> setaccio;
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n
    mpz_class base = sqrt(n) + 1;
    initializeSieve(n, base, L, factorBase, setaccio);

    #ifdef DEBUG
    printSetaccio(setaccio, base);
    #endif



    // controllare che il setaccio sia stato inizializzato correttamente
    // fatto con un test



    // attivare sieve

        // durante il setaccio i primi/le potenze vengono aggiunti a un numero dal più grande al più piccolo (mentre si analizzano i numeri precedenti)
        // --> se si inseriscono con push_front poi si ritrovano dal più piccolo al più grande (mentre si analizza quel numero)

    std::vector<smoothElem> smooths;
    unsigned long long missing = 0LL + numPrimes + mpz_sizeinbase(n.get_mpz_t(), 2), numEqs;
    mpz_class a = -1;   // dato che nella funzione l'incremento avviene prima del resto, in questo modo il primo valore di a sarà 0
    mpz_class baseSquaredMinusN = base*base - n, twoBase = base<<1;

    numEqs = activateSieve(missing, logMaxP2, a, base, baseSquaredMinusN, twoBase, factorBase, setaccio, smooths);

    #ifdef DEBUG
    printSmooths(smooths);
    #endif

    // preprocessing: calcolare exponent vectors, fare (o no) filtering, se non ci sono abbastanza smooth numbers tornare al sieving
        // dato ogni numero in smooths:
            // tirare via le potenze dei primi noti segnandosi quali sono quelli con esp dispari (in un qualche ordine)
            // segnarsi ciò che rimane (se diverso da 1) come nuovo primo che compare con esp 1 (quindi segnarlo nell'exp vector ridotto) e aggiungerlo ad un set di nuovi primi
        // se il numero di numeri trovati non supera numPrimes + numNewPrimes riprendere sieving, altrimenti si può ripulire il setaccio, ma segnandosi l'ultimo elem visto come nuova base
        // fare un unico set con tutti i primi che hanno degli exp dispari, e per ognuno di questi ricavare la posizione ordinata,
            // così da poter fare poi gli exponent vectors con le posizioni dei primi (notazione compatta per matrice densa)




    // algebra lineare





    // differenza quadrati: radici e gcd


    // in caso di fail...

}




unsigned long long activateSieve(unsigned long long toFind, unsigned short maxLogp2, mpz_class& a, mpz_class& base, mpz_class& baseSquaredMinusN, mpz_class& twoBase,  std::unordered_map<mpz_class, unsigned short>& factorBase, std::unordered_map<mpz_class, elemSetaccio>& setaccio, std::vector<smoothElem>& smooths){   // il valore di ritorno è il numero di smooth trovati
    mpz_class x, y, P;
    unsigned long long found = 0;
    unsigned short lastLog = 0;
    elemSetaccio divisors;

    while(found < toFind){
        a++;                    // a era l'ultimo visto

        auto it = setaccio.find(a);
        if(it == setaccio.end()) continue;

        // ora it->second dovrebbe essere la lista con potenze dei primi che dividono y(x(a)) e relativi logaritmi
        divisors.swap(it->second);
        setaccio.erase(it);

        // ora divisors dovrebbe essere la lista con potenze dei primi che dividono y(x(a)) e relativi logaritmi

        unsigned short sum_of_logs = 0, l = 0;
        for(auto& coppia : divisors){
            P = coppia.first;
            l = coppia.second;

            // passare avanti la potenza che divide y
            auto j = setaccio.find(a + P);
            if(j == setaccio.end()){
                setaccio[a+P] = elemSetaccio{std::make_pair(P, l)};
            } else{
                j->second.push_front(std::make_pair(P, l));
            }

            // sommare i log
            sum_of_logs += l;
        }

        if(sum_of_logs < lastLog){
            divisors.clear();
            continue;
        }

        y = baseSquaredMinusN + a*(twoBase + a);
        double logyD = log2(y.get_d());
        unsigned short logy = (unsigned short) (logyD * (1<<6));
        // lastLog is < logy to allow inclusion of "large" prime    - da Wikiversity

        lastLog = logy - maxLogp2;
        if(sum_of_logs < lastLog){
            divisors.clear();
            continue;
        }

        // abbiamo un elemento probabilmente smooth, tranne per al più un large prime
        found++;
        smooths.push_back(smoothElem(a+base, y));
        smoothElem & elem = smooths.back();

        while(!divisors.empty()){
            P = (divisors.front()).first;
            divisors.pop_front();

            if(factorBase.count(P)){    // userei contains ma c'è solo da c++20
                elem.primes.push_front(P);
            }
        }
    }

    return found;
}




void initializeSieve(const mpz_class& n, mpz_class& base, mpz_class& L, std::unordered_map<mpz_class, unsigned short>& factorBase, std::unordered_map<mpz_class, elemSetaccio>& setaccio){
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n

    mpz_class p, Pot, Pot_, a, tmp1, tmp2, tmp3;
    std::vector<mpz_class> roots;
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

        insertRoots(roots, base, p, setaccio, a, l);


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

            insertRoots(roots, base, Pot, setaccio, a, l);

            Pot_ = Pot;
            Pot *= p;
        }
    }
}



void insertRoots(std::vector<mpz_class>& roots, mpz_class& base, mpz_class& P, std::unordered_map<mpz_class, elemSetaccio>& setaccio, mpz_class& a, unsigned short l){
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
        }
    }
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
    return logMaxP << 1;
}



void choose_params(mpz_class &n, mpz_class &B, mpz_class &L){
    // B = exp((1/2 + o(1))*(logn loglogn)^(1/2)), e controllare che sia almeno 5
        // (altrimenti ci vuole troppo tempo per i numeri piccoli perché non si trovano abbastanza B-smooth)

    // L = exp((logn loglogn)^(1/2)-(1/2)*loglogn), e controllare che sia almeno 100

    double logn = log(n.get_d()), loglogn = log(logn);
    double  b = exp((1/2.0 + 3/logn)*sqrt(logn * loglogn)),                     // scelto 3/logn come o(1)
            l = exp((1.0 + 6/logn)*sqrt(logn * loglogn) - (1/2.0)*loglogn);     // qui usato il doppio perché c'è circa un fattore 2 di differenza

    B = b; L = l;
    if(B < 5) B=5;
    if(L < 100) L=100;
}
