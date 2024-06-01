#include "quadratic_sieve.hpp"
//#include<iostream>

void quadratic_sieve(mpz_class& n, mpz_class& p, mpz_class& q){

    // scelta parametri
    // B e magari approx lunghezza intervallo di sieving
    mpz_class B, L;
    //std::tie(B, L) = choose_params(n);
    choose_params(n, B, L);
    //std::cout<<B<<" "<<L<<std::endl;

    // costruire factor base



    // inizializzare sieve



    // attivare sieve



    // preprocessing: calcolare exponent vectors, fare (o no) filtering, se non ci sono abbastanza smooth numbers tornare al sieving



    // algebra lineare



    // differenza quadrati: radici e gcd


}




void choose_params(mpz_class &n, mpz_class &B, mpz_class &L){
    // B = exp((1/2)*(logn loglogn)^(1/2)), e controllare che sia almeno 5
        // (altrimenti ci vuole troppo tempo per i numeri piccoli perché non si trovano abbastanza B-smooth)

    // L = exp((logn loglogn)^(1/2)-(1/2)*loglogn), e controllare che sia almeno 100

    double logn = log(n.get_d()), loglogn = log(logn);
    double b = exp((1/2.0)*sqrt(logn * loglogn)), l = exp(sqrt(logn * loglogn) - (1/2.0)*loglogn);

    B = b; L = l;
    //return std::make_pair(B, L);
}