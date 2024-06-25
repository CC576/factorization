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
#include <sys/resource.h>

#include "src/utils/statsVar.hpp"
#include "src/utils/utils_mpz/hash_mpz.hpp"
#include "src/utils/utils_mpz/mpz_ull.hpp"
#include "src/utils/utils_modP/roots_modP.hpp"

#include "src/trial_div/trial_division.hpp"
#include "src/fermat/fermat.hpp"
#include "src/quadratic_sieve/quadratic_sieve.hpp"
#include "src/gnfs/gnfs.hpp"

//#include <NTL/ZZ.h>
//#include <NTL/ZZ_pXFactoring.h>

/*
#include <linbox/ring/modular.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/algorithms/mg-block-lanczos.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/vector/vector.h>
using Field = Givaro::UnparametricZRing<int>;
*/

//void testGmp();

int main(int argc, char** argv){
    // l'algoritmo da usare dev'essere argv[1], n dev'essere argv[2]

    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " <algorithm> <n>";
        return 1;
    }

    mpz_class n;
    n = argv[2];
    mpz_class p,q;
    //std::cout<<mpz_scan1(n.get_mpz_t(), 0)<<std::endl;
    int alg = argv[1][0] - '0';

    if(argc > 3) printStats = true;

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

    case 4:
        gnfs(n, p, q);
        break;

    default:
        std::cerr << "Invalid algorithm " << alg << std::endl;
        return 2;
    }

    std::cout << "(" << p << ", " << q << ")" << std::endl;
    struct rusage usage;
    int ret;
    ret = getrusage(RUSAGE_SELF, &usage);
    std::cout   << "\"userTime\": [" << usage.ru_utime.tv_sec << ", " << usage.ru_utime.tv_usec << "], "
                << "\"systemTime\": [" << usage.ru_stime.tv_sec << ", " << usage.ru_stime.tv_usec << "], "
                << "\"RSS (in kB)\": " << usage.ru_maxrss << std::endl;

    return 0;
}



/*
#include "src/gnfs/factorBases.hpp"
#include "src/utils/utils_mpz/mpz_ZZ.hpp"
#include "src/utils/utils_NTL/utils_ZZX.hpp"
#include <NTL/ZZ_pXFactoring.h>
int main(int argc, char** argv){
    mpz_class nMPZ; // = 1497130049;    // 3221763169
    nMPZ = argv[2];
    ZZ n;
    mpz_2_ZZ(nMPZ, n);
    long d;
    ZZ m, B;
    ZZX f;
    chooseParams(nMPZ, d, m, f, B);

    std::cout << "f: " << f << std::endl;

    std::vector<std::pair<long, uint8_t>> primes;
    genPrimesList(primes, B);

    ZZ_pX fp;

    for(auto& [p, l] : primes){
        ZZ_p::init(conv<ZZ>(p));
        fp = conv<ZZ_pX>(f);

        Vec<Pair<ZZ_pX, long>> factors;
        CanZass(factors, fp);

        for(auto& [fattore, esp] : factors){
            if(deg(fattore) == 1 && esp > 1){
                std::cout << "Found multiple root: " << fattore << " with multiplicity " << esp << " mod " << p << std::endl;
                MakeMonic(fattore);
                ZZ_p coso = -coeff(fattore, 0);
                ZZ coso2 = conv<ZZ>(coso);
                ZZ tmp; ZZX_eval(tmp, f, coso2);
                std::cout << tmp << ", e: ";

                uint8_t e = (uint8_t) 0;
                if(tmp == 0){
                    e = (uint8_t) 255;
                    tmp = 1;
                }
                while(tmp%p == 0){
                    e++; tmp/=p;
                }
                std::cout <<  int(e) << std::endl;
            }
        }
        //std::cout << factors << std::endl;
    }

}//*/




/*int main6(){          // printa anche p=2
    PrimeSeq s;
   long p;

   p = s.next();
   while (p <= 1000) {
      std::cout << p << "\n";
      p = s.next();
   }
}*/




/*int main5(){

    NTL::ZZ p;
    //std::cin >> p;
    p = 1531;
    NTL::ZZ_p::init(p);

    NTL::ZZ_pX f(2, 1);
    NTL::SetCoeff(f, 1, -4);
    NTL::SetCoeff(f, 0, 3);
    //f.ZZ_pX::ZZ_pX.se
    //std::cin >> f;

    NTL::Vec< NTL::Pair< NTL::ZZ_pX, long > > factors;

    NTL::CanZass(factors, f);  // calls "Cantor-Zassenhaus" algorithm

    std::cout << factors << "\n";

    return 0;
}*/


/*int main4(){

    Field f;
    LinBox::MGBlockLanczosSolver<Field> solver(f, LinBox::Method::BlockLanczos());

    LinBox::SparseMatrix<Field> mat(f, 5, 5);
    mat.setEntry(0, 0, 1); // Set entry at (0, 0) to 1
    mat.setEntry(1, 2, 3); // Set entry at (1, 2) to 3
    mat.setEntry(4, 4, 6); // Set entry at (4, 4) to 6

    LinBox::DenseMatrix<Field> nullspaceMatrix(f, 5, 5);

    int solcount = solver.sampleNullspace(mat, nullspaceMatrix);
    std::cout << solcount << std::endl;
    std::cout << nullspaceMatrix << std::endl;

    return 0;
}*/

int main3(){
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


/*int main2(){
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

    return 0;
}*/