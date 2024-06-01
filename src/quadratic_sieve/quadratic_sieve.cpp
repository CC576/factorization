#include "quadratic_sieve.hpp"
#include<iostream>

void quadratic_sieve(mpz_class& n, mpz_class& fattore1, mpz_class& fattore2){

    // scelta parametri
    // B e magari approx lunghezza intervallo di sieving
    mpz_class B, L;
    choose_params(n, B, L);
    //std::cerr << B << std::endl << std::endl;

    // costruire factor base
    std::vector<std::pair<mpz_class, unsigned short>> factorBase;
    unsigned short logMaxP2 = buildFactorBase(n, B, factorBase);
    unsigned long long numPrimes = factorBase.size();

    /*std::cerr << std::endl << std::endl << factorBase.size() << std::endl << std::endl;
    for(auto coppia : factorBase){
        std::cerr << coppia.first << std::endl;
    }*/


    // inizializzare sieve
    std::unordered_map<mpz_class, elemSetaccio> setaccio;
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n
    mpz_class base = sqrt(n);
    initializeSieve(n, base, L, factorBase, setaccio);



    // controllare che il setaccio sia stato inizializzato correttamente




    // attivare sieve





    // preprocessing: calcolare exponent vectors, fare (o no) filtering, se non ci sono abbastanza smooth numbers tornare al sieving





    // algebra lineare





    // differenza quadrati: radici e gcd


}


unsigned long getOrdInSylow(const mpz_class& p, mpz_class& u, unsigned long s, mpz_class& tmp){
    unsigned long rk = 0;   // esponente da dare a 2 per avere u^(2^rk) == 1 % p
    tmp = u;
    while(tmp != 1 && rk < s){
        tmp = (tmp*tmp)%p;
        rk++;
    }
    return rk;
}



void findRoot(mpz_class& a, const mpz_class& p, mpz_class& v, mpz_class& u, mpz_class& m, mpz_class& c){    // andrebbero bene anche dei long long probabilmente
    // implementa Tonelli-Shanks
std::cerr<<"here"<<std::endl<<std::flush;
    if(p%4 == 3){
        m = (p+1)/4;
        mpz_powm(v.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
        return;
    }
std::cerr<<"here2\n";
    m = p-1;
    unsigned long s = mpz_scan1(m.get_mpz_t(), 0), rk;  // s è l'esponente di 2 più grande in p-1
std::cerr<<s<<"\n"<<std::flush;
    m >>= s;     // adesso m è la parte dispari di p-1
std::cerr<<s<<" "<<m<<"\n"<<std::flush;
    mpz_powm(u.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());   // u == a^m mod p
    m = (m+1)/2;
    mpz_powm(v.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());   // v == a^((m+1)/2) mod p
    m = 2*m - 1;

    if(u==1) return;

    getNonQResidue(p, c);
    mpz_powm(c.get_mpz_t(), c.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());   // da qui m non dovrebbe servire più; c è nel sylow subgroup

    while(u != 1){
        rk = getOrdInSylow(p, u, s, m);
        m = 2;
        mpz_pow_ui(m.get_mpz_t(), m.get_mpz_t(), s-rk-1);   // m = 2^(s-rk-1)

        mpz_powm(c.get_mpz_t(), c.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
        v = (v*c) % p;
        c = (c*c) % p;
        u = (u*c) % p;

        s = rk; // i nuovi elementi hanno ord che divide 2^rk
    }
}



void initializeSieve(mpz_class& n, mpz_class& base, mpz_class& L, std::vector<std::pair<mpz_class, unsigned short>>& factorBase, std::unordered_map<mpz_class, elemSetaccio>& setaccio){
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

        for(auto& r : roots){
            // individuale il valore positivo di a più piccolo per cui a+base == r mod p
            // --> a == r-base mod p
            a = (r - base) % p; // la divisione tronca verso 0
            if(a < 0) a += p;

            auto it = setaccio.find(a);
            if(it == setaccio.end()){
                setaccio[a] = elemSetaccio{std::make_pair(p, l)};
            }else{
                it->second.push_front(std::make_pair(p, l));
            }
        }


        //potenze (ricordarsi che per p=2 c'è un'unica radice, per Pot = 4 ce ne sono 2 per quadrato, da Pot=8 in poi ogni quadrato ha  4 radici)


    }
}


void getNonQResidue(const mpz_class& p, mpz_class& c){
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    int l;
    do{
        mpz_urandomm (c.get_mpz_t(), state, p.get_mpz_t());
        l = mpz_legendre(c.get_mpz_t(), p.get_mpz_t());
    } while(l != -1);
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

        if(p != 2 && mpz_legendre(n.get_mpz_t(), tmp.get_mpz_t()) != 1) continue;

        double log2pD = log2(double(p));
        unsigned short log2p = (unsigned short) (log2pD * (1<<6));

        factorBase.push_back(std::make_pair(tmp, log2p));

        logMaxP = std::max(logMaxP, log2p);
    }
    return logMaxP << 1;
}



void choose_params(mpz_class &n, mpz_class &B, mpz_class &L){
    // B = exp((1/2)*(logn loglogn)^(1/2)), e controllare che sia almeno 5
        // (altrimenti ci vuole troppo tempo per i numeri piccoli perché non si trovano abbastanza B-smooth)

    // L = exp((logn loglogn)^(1/2)-(1/2)*loglogn), e controllare che sia almeno 100

    double logn = log(n.get_d()), loglogn = log(logn);
    double b = exp((1/2.0)*sqrt(logn * loglogn)), l = exp(sqrt(logn * loglogn) - (1/2.0)*loglogn);

    B = b; L = l;
    if(B < 5) B=5;
    if(L < 100) L=100;
}