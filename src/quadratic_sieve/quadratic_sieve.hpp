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
#include "../utils/utils_modP/roots_modP.hpp"

typedef std::forward_list<std::pair<mpz_class, unsigned short>> elemSetaccio;

void quadratic_sieve(mpz_class&, mpz_class&, mpz_class&);

void choose_params(mpz_class&, mpz_class&, mpz_class&);
unsigned short buildFactorBase(mpz_class&, mpz_class&, std::unordered_map<mpz_class, unsigned short>&);
void initializeSieve(const mpz_class&, mpz_class&, mpz_class&, std::unordered_map<mpz_class, unsigned short>&, std::unordered_map<mpz_class, elemSetaccio>&);
void insertRoots(std::vector<mpz_class>&, mpz_class&, mpz_class&, std::unordered_map<mpz_class, elemSetaccio>&, mpz_class& a, unsigned short);



#endif