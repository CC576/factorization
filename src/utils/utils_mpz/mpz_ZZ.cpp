#include "mpz_ZZ.hpp"
#include <string>
#include <sstream>
//#include <iostream>


void mpz_2_ZZ(const mpz_t& a, ZZ& B){
    std::string s;
    s.resize(mpz_sizeinbase(a, 10) + 2);
    mpz_get_str(s.data(), 10, a);
    B = conv<ZZ>(s.c_str());
}

void mpz_2_ZZ(const mpz_class& A, ZZ& B){
    std::string s;
    s = A.get_str();
    B = conv<ZZ>(s.c_str());
}



void mpz_from_ZZ(const ZZ& A, mpz_t& b){
    mpz_init(b);

    std::stringstream coso;
    coso << A;

    mpz_set_str(b, coso.str().c_str(),10);
}

void mpz_from_ZZ(const ZZ& A, mpz_class& B){
    std::stringstream coso;
    coso << A;
    coso >> B;
}