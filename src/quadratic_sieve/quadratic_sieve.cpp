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
    std::vector<std::pair<mpz_class, unsigned short>> factorBase;
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




    // attivare sieve





    // preprocessing: calcolare exponent vectors, fare (o no) filtering, se non ci sono abbastanza smooth numbers tornare al sieving





    // algebra lineare





    // differenza quadrati: radici e gcd


}




void initializeSieve(const mpz_class& n, mpz_class& base, mpz_class& L, std::vector<std::pair<mpz_class, unsigned short>>& factorBase, std::unordered_map<mpz_class, elemSetaccio>& setaccio){
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n

    mpz_class p, Pot, Pot_, a, tmp1, tmp2, tmp3;
    std::vector<mpz_class> roots;
    for(auto& coppia : factorBase){
        p = coppia.first;
        unsigned short l = coppia.second;

        if(p == 2){
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



unsigned short buildFactorBase(mpz_class& n, mpz_class& B, std::vector<std::pair<mpz_class, unsigned short>>& factorBase){
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

        factorBase.push_back(std::make_pair(tmp, log2p));

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
