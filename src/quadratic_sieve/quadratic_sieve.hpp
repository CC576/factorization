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
#include <unordered_map>
#include <forward_list>
#include "../utils/utils_mpz/hash_mpz.hpp"
#include "../utils/utils_mpz/mpz_ull.hpp"

void quadratic_sieve(mpz_class&, mpz_class&, mpz_class&);

void choose_params(mpz_class&, mpz_class&, mpz_class&);
void buildFactorBase(mpz_class&, mpz_class&, std::vector<std::pair<mpz_class, unsigned short>>&);
void initializeSieve(mpz_class&, mpz_class&, mpz_class&, std::vector<std::pair<mpz_class, unsigned short>>&, std::unordered_map<unsigned int, std::forward_list<std::pair<mpz_class, unsigned short>>>&);
unsigned int findRoot(unsigned int, unsigned int);


#endif