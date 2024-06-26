#include "hash_ZZ.hpp"

namespace std {

size_t hash<NTL::ZZ>::operator()(const NTL::ZZ& seed) const {
    long nb = NTL::NumBytes(seed);

    string buf;
    buf.resize(nb);
    NTL::BytesFromZZ(reinterpret_cast<unsigned char*>(&buf[0]), seed, nb);

    size_t result = string_hasher(buf);
    return result;
}

size_t hash<pair<NTL::ZZ, NTL::ZZ>>::operator()(const pair<NTL::ZZ, NTL::ZZ>& seed) const {
    size_t h1 = zz_hasher(seed.first);
    size_t h2 = zz_hasher(seed.second);
    return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}

}

