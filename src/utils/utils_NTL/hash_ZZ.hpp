#ifndef HASH_ZZ_HPP_
#define HASH_ZZ_HPP_

#include <cstddef>
#include <string>
#include <NTL/ZZ.h>

namespace std {

template <>
struct hash<NTL::ZZ> {
    hash<string> string_hasher{};
    size_t operator()(const NTL::ZZ& seed) const;
};

template <>
struct hash<pair<NTL::ZZ, NTL::ZZ>> {
    hash<NTL::ZZ> zz_hasher{};
    size_t operator()(const pair<NTL::ZZ, NTL::ZZ>& seed) const;
};

}

#endif /* HASH_MPZ_HPP_ */

