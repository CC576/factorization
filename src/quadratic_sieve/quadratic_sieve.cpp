#include "quadratic_sieve.hpp"
//#include<iostream>

void quadratic_sieve(mpz_class& n, mpz_class& fattore1, mpz_class& fattore2){

    // scelta parametri
    // B e magari approx lunghezza intervallo di sieving
    mpz_class B, L;
    choose_params(n, B, L);
    //std::cerr << B << std::endl << std::endl;

    // costruire factor base
    std::vector<std::pair<mpz_class, unsigned short>> factorBase;
    buildFactorBase(n, B, factorBase);
    unsigned long long numPrimes = factorBase.size();

    /*std::cerr << std::endl << std::endl << factorBase.size() << std::endl << std::endl;
    for(auto coppia : factorBase){
        std::cerr << coppia.first << std::endl;
    }*/


    // inizializzare sieve
    //std::unordered_map<unsigned long long, std::forward_list<std::pair<mpz_class, unsigned short>>> setaccio;
    // il setaccio sarà indicizzato da a=1..., con a=x-floor(sqrt(n)) --> y(x)=(a+base)^2 - n
    mpz_class base = sqrt(n);
    //initializeSieve(n, base, L, factorBase, setaccio);



    // controllare che il setaccio sia stato inizializzato correttamente




    // attivare sieve





    // preprocessing: calcolare exponent vectors, fare (o no) filtering, se non ci sono abbastanza smooth numbers tornare al sieving





    // algebra lineare





    // differenza quadrati: radici e gcd


}


void findRoot(mpz_class& n, mpz_class& p, mpz_class& r){

}




void initializeSieve(mpz_class& n, mpz_class& base, mpz_class& L, std::vector<std::pair<mpz_class, unsigned short>>& factorBase, std::unordered_map<unsigned int, std::forward_list<std::pair<mpz_class, unsigned short>>>& setaccio){
    mpz_class p, Pot;
    std::vector<mpz_class> roots;
    for(auto& coppia : factorBase){
        if(p == 2){
            roots.push_back(0);
        } else{
            roots.resize(2);
            findRoot(n, p, roots[0]);
            roots[1] = p-roots[0];
        }

        for(auto& r : roots){

        }


        //potenze (ricordarsi che per p=2 c'è un'unica radice)
    }
}






void buildFactorBase(mpz_class& n, mpz_class& B, std::vector<std::pair<mpz_class, unsigned short>>& factorBase){
    // effettua anche controllo simbolo di Legendre, che per versione MPQS va spostato nell'inizializzazione del sieve

    mpz_class tmp;
    unsigned long long b = mpz_2_ull(B, tmp);  // sto assumendo che B sia < 2^64 e che possa stare in un long long

    std::vector<bool> composite(b+1);

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
    }
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