#ifndef QUADRATIC_HPP
#define QUADRATIC_HPP

#include <gmp.h>
#include <gmpxx.h>
//#include <mpfr.h>
#include <cmath>
#include <blanczos.h>
#include <utility>
#include <tuple>
#include <vector>
#include <forward_list>
#include <unordered_map>
//#include <unordered_set>
//#include <bit>
#include <bitset>
#include "../utils/utils_mpz/hash_mpz.hpp"
#include "../utils/utils_mpz/mpz_ull.hpp"
#include "../utils/utils_modP/roots_modP.hpp"

typedef std::forward_list<std::pair<mpz_class, unsigned short>> elemSetaccio;   // scelta infelice per questo nome, sarebbe più corretto valueSetaccio
//typedef std::vector<std::pair<mpz_class, unsigned short>> elemSetaccio;
struct smoothElem{
  mpz_class x;
  mpz_class y;
  elemSetaccio primes;

  smoothElem(mpz_class x_, mpz_class y_){
    x = x_; y = y_;
    primes = std::forward_list<std::pair<mpz_class, unsigned short>>();
    //primes = std::vector<std::pair<mpz_class, unsigned short>>(0);
  }
};


void quadratic_sieve(mpz_class&, mpz_class&, mpz_class&);

void choose_params(mpz_class&, mpz_class&, mpz_class&);
unsigned short buildFactorBase(mpz_class&, mpz_class&, std::unordered_map<mpz_class, unsigned short>&);
void initializeSieve(const mpz_class&, mpz_class&, mpz_class&, std::unordered_map<mpz_class, unsigned short>&, std::unordered_map<mpz_class, elemSetaccio>&);
void insertRoots(std::vector<mpz_class>&, mpz_class&, mpz_class&, std::unordered_map<mpz_class, elemSetaccio>&, mpz_class& a, unsigned short);
unsigned long long activateSieve(unsigned long long, unsigned short, mpz_class&, mpz_class&, mpz_class&, mpz_class&, /*std::unordered_map<mpz_class, unsigned short>&,*/ std::unordered_map<mpz_class, elemSetaccio>&, std::vector<smoothElem>&, unsigned short&);
unsigned long long getMatrix(std::vector<smoothElem>&, std::vector<uint32_t>&, std::unordered_map<mpz_class, uint32_t>&);

#endif