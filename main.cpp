#include <stdio.h>
#include <iostream>
#include <blanczos.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
//#include <Singular/libsingular.h>
//#include <Polynomial/Polynomial.hpp>
#include <tuple>
#include <cassert>

#include "src/utils/utils_mpz/hash_mpz.hpp"
#include "src/utils/utils_mpz/mpz_ull.hpp"
#include "src/utils/utils_modP/roots_modP.hpp"

#include "src/trial_div/trial_division.hpp"
#include "src/fermat/fermat.hpp"
#include "src/quadratic_sieve/quadratic_sieve.hpp"

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
    //std::cout<<mpz_scan1(n.get_mpz_t(), 0)<<std::endl;
    int alg = argv[1][0] - '0';

    switch (alg)
    {
    case 1:
        trial_division(n, p, q);
        break;

    case 2:
        fermat(n, p, q);
        break;

    case 3:
        quadratic_sieve(n, p, q);
        break;

    default:
        std::cerr << "Invalid algorithm " << alg << std::endl;
        return -2;
    }

    std::cout << "(" << p << ", " << q << ")" << std::endl;

    return 0;
}

int main2(){
    /*mpz_class coso, coso2;
    coso = "18446744073709551616";
    coso2 = coso>>32;
    std::cerr<<coso2<<std::endl;
    return 0;*/

    mpz_class root = 102, p = 257, a = (root*root)%p, tmp1, tmp2, tmp3, v;
    findRoot(a, p, v, tmp1, tmp2, tmp3);
    assert((v*v)%p == a);

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