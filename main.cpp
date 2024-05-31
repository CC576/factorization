#include <stdio.h>
#include <iostream>
#include <blanczos.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
//#include <Singular/libsingular.h>
//#include <Polynomial/Polynomial.hpp>
#include <tuple>

#include "src/trial_div/trial_division.hpp"

//void testGmp();

int main(int argc, char** argv){
    // l'algoritmo da usare dev'essere argv[1], n dev'essere argv[2]

    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " <algorithm> <n>";
        return -1;
    }

    mpz_class n;
    n = argv[2];

    mpz_class p,q;

    std::tie(p, q) = trial_division(n);

    std::cout << "(" << p << ", " << q << ")" << std::endl;

    return 0;
}



/*void testGmp(){
    //std::float16_t A;
    uint32_t* B;
    uint64_t N=0;
    uint32_t Nrow=1, Ncol=65;
    uint64_t* result;
    //blanczos(B, N, Nrow, Ncol, result);
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
}*/