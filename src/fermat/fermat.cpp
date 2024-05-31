#include "fermat.hpp"
#include <iostream>


std::pair<mpz_class, mpz_class> fermat(mpz_class n){
    //std::cerr << n << std::endl;
    mpz_class base = sqrt(n), h = base*base - n, i = 0;
    mpz_class x, y = 0, nHalves = n/2, twoBase = 2*base;

    while(y < nHalves){
        h = h + twoBase + 2*i + 1;
        i = i+1;
        y = sqrt(h);

        if(y*y == h){
            x = base + i;
            break;
        }
    }

    mpz_class p = x-y, q = x+y;

    return std::make_pair(p, q);
}