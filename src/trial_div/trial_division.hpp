#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <gmp.h>
#include <gmpxx.h>
#include <utility>

std::pair<mpz_class, mpz_class> trial_division(mpz_class);

#endif