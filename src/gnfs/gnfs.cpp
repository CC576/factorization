#include "gnfs.hpp"

#ifdef DEBUG
#include "../utils/utils_print/print_stuff.hpp"
#include <cassert>
#include<iostream>
#endif
#include<iostream>
#include<sys/resource.h>

#include "../utils/utils_mpz/mpz_ZZ.hpp"



void gnfs(const mpz_class& nMPZ, mpz_class& fattore1, mpz_class& fattore2){


    // 1. preparazione

    // 1.a convertire n in uno ZZ
    ZZ n;
    mpz_2_ZZ(nMPZ, n);
    std::cout << n << std::endl;

    /*mpz_class coso;
    mpz_from_ZZ(n, coso);
    std::cout << coso << std::endl;*/   // la conversione sembra funzionare

    // 1.b scelta parametri: costruire polinomio (e testare irriducibilità)


    // 1.c costruire factor bases




    // 2. sieving

    // 2.a sieving razionale

    // 2.b sieving algebrico




    // 3. dipendenze lineari

    // 3.a preprocessing (trial division, segno, e quadratic characters)

    // 3.b block lanczos




    // iterare sulle dipendenze lineari trovate

    // 4. radici

    // 4.a radice razionale (includendo f'(m))

    // 4.b trovare finite fields applicabili (volendo si può fare durante costruzione factor bases)/aggiungerne altri se non bastano

    // 4.c radice in finite fields, metodo di Couveignes

    // 4.d CRT




    // 5. fattorizzare n

    // 5.a GCD

    // 5.b convertire p e q in mpz_class per metterli in fattore1 e fattore2




    // 5. gestire fail...
}