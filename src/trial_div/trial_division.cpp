#include "trial_division.hpp"


void trial_division(mpz_class& n, mpz_class& p, mpz_class& q){
    mpz_class base = sqrt(n);

    for(p = 2; p <= base; p++){
        if(n%p == 0){
            q = n/p;
            break;
        }
    }
}