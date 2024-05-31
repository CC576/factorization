#include "trial_division.hpp"


std::pair<mpz_class, mpz_class> trial_division(mpz_class n){
    mpz_class base = sqrt(n);

    mpz_class p,q;
    for(p = 2; p <= base; p++){
        if(n%p == 0){
            q = n/p;
            break;
        }
    }

    return std::make_pair(p, q);
}