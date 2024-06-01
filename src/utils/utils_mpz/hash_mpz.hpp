#ifndef HASH_MPZ_HPP_
#define HASH_MPZ_HPP_

#include <gmpxx.h>

namespace std {

template<> struct hash<mpz_srcptr> {
    size_t operator()(const mpz_srcptr x) const;
};

template<> struct hash<mpz_t> {
    size_t operator()(const mpz_t x) const;
};

template<> struct hash<mpz_class> {
    size_t operator()(const mpz_class &x) const;
};

}

#endif /* HASH_MPZ_HPP_ */
