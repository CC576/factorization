#include <stdio.h>
#include <iostream>
#include <blanczos.h>
#include <gmp.h>
#include <gmpxx.h>
//#include <Singular/libsingular.h>
//#include <Polynomial/Polynomial.hpp>

int main(int, char**){
    //std::float16_t A;
    uint32_t* B;
    uint64_t N=0;
    uint32_t Nrow=1, Ncol=65;
    uint64_t* result;
    blanczos(B, N, Nrow, Ncol, result);
    printf("Hello, from factoring_algorithms!\n");

    mpz_class v;
    v = "-42384675928734659287364059287340958720938475029384750293845";
    std::cout << v << std::endl;
    mpz_class w = v*v;
    mpz_class x = v-1;
    std::cout << w << std::endl;

    mpz_class y;
    mpz_powm(x.get_mpz_t(), v.get_mpz_t(), w.get_mpz_t(), x.get_mpz_t());
    std::cout << y << std::endl;
}
