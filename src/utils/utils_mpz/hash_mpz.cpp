#include "hash_mpz.hpp"
#include <cstddef>
#include <string_view>

constexpr size_t pi_size_t() {
    if (sizeof(size_t) == 4) {
        return 0xc90fdaa2; // floor(pi/4 * 2^32)
    } else if (sizeof(size_t) == 8) {
        return 0xc90fdaa22168c234; // floor(pi/4 * 2^64)
    } else {
        throw std::logic_error(
                "sizeof(size_t) not supported. only 32 or 64 bits are supported. you can easily add the required code for other sizes.");
    }
}

inline size_t scramble(size_t v) {
    return v ^ (pi_size_t() + (v << 6) + (v >> 2));
}

namespace std {

size_t std::hash<mpz_srcptr>::operator()(const mpz_srcptr x) const {
    string_view view { reinterpret_cast<char*>(x->_mp_d), abs(x->_mp_size)
            * sizeof(mp_limb_t) };
    size_t result = hash<string_view> { }(view);

    // produce different hashes for negative x
    if (x->_mp_size < 0) {
        result = scramble(result);
    }

    return result;
}

size_t hash<mpz_t>::operator()(const mpz_t x) const {
    return hash<mpz_srcptr> { }(static_cast<mpz_srcptr>(x));
}

size_t hash<mpz_class>::operator()(const mpz_class &x) const {
    return hash<mpz_srcptr> { }(x.get_mpz_t());
}

}