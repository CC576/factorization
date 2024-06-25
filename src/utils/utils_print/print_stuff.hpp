#ifndef PRINT_STUFF_HPP
#define PRINT_STUFF_HPP

#include<iostream>

#include <gmp.h>
#include <gmpxx.h>

#include <utility>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <forward_list>
//#include "../utils_mpz/hash_mpz.hpp"
#include "../../quadratic_sieve/quadratic_sieve.hpp"
#include "../../gnfs/gnfs.hpp"
/*
typedef std::forward_list<std::pair<mpz_class, unsigned short>> elemSetaccio;
typedef struct smoothElem_{       // rischio di fare confusione con queste definizioni di tipi sparse in più file
  mpz_class x;
  mpz_class y;
  std::forward_list<mpz_class> primes;

  smoothElem_(mpz_class x_, mpz_class y_){
    x = x_; y = y_;
    primes = std::forward_list<mpz_class>();
  }
} smoothElem__;
*/
void printFactorBase(std::unordered_map<mpz_class, unsigned short>&);
void printSetaccio(std::unordered_map<mpz_class, elemSetaccio>&, mpz_class&);
void printSmooths(std::vector<smoothElem> &);
void printUsedPrimes(std::unordered_map<mpz_class, uint32_t> &);

void printFBgnfs(factorBase& FB);

#endif