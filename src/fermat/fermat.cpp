#include "fermat.hpp"
//#include <iostream>


void fermat(mpz_class& n, mpz_class& p, mpz_class& q){
    //std::cerr << n << std::endl;
    mpz_class base = sqrt(n), h = base*base - n;
    //mpz_class i = 0;
    mpz_class x = base, y = 0, nHalves = n/2;
    //mpz_class twoBase = 2*base;

    while(y < nHalves){
        //h = h + twoBase + 2*i + 1;
        //i = i+1;
        h = h + (x<<1) + 1;
        x++;
        y = sqrt(h);

        if(y*y == h){
            //x = base + i;
            break;
        }
    }

    mpz_class p = x-y, q = x+y;
}