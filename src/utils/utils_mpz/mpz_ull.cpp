#include "mpz_ull.hpp"


unsigned long long mpz_2_ull(mpz_t A, mpz_t B){

    unsigned long long n = mpz_get_ui(A) & ((1LL << 32)-1);
    mpz_mul_2exp(B, A, -32);
    n += ((unsigned long long) mpz_get_ui(B)) << 32;

    return n;
}

unsigned long long mpz_2_ull(mpz_class& A, mpz_class& B){
    unsigned long long n = A.get_ui() & ((1LL << 32)-1);
    B = A >> 32;
    n += ((unsigned long long) B.get_ui()) >> 32;

    return n;
}


void mpz_from_ull(unsigned long long n, mpz_t A){
    mpz_set_ui(A, (unsigned int) (n >> 32));
    mpz_mul_2exp(A, A, 32);
    mpz_add_ui(A, A, (unsigned int) (n & ((1LL << 32)-1)));
}
void mpz_from_ull(unsigned long long n, mpz_class& A){
    A = (unsigned int) (n >> 32);
    A <<= 32;
    A += (unsigned int) (n & ((1LL << 32)-1));
}