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
#include <forward_list>
#include "../utils_mpz/hash_mpz.hpp"

typedef std::forward_list<std::pair<mpz_class, unsigned short>> elemSetaccio;

void printFactorBase(std::vector<std::pair<mpz_class, unsigned short>>&);
void printSetaccio(std::unordered_map<mpz_class, elemSetaccio>&, mpz_class&);

#endif