#include "roots_modP.hpp"




void findRoot(const mpz_class& a, const mpz_class& p, mpz_class& v, mpz_class& u, mpz_class& m, mpz_class& c){    // andrebbero bene anche dei long long probabilmente
    // implementa Tonelli-Shanks
    //std::cerr<<"here"<<std::endl<<std::flush;
    if(p%4 == 3){
        m = (p+1)/4;
        mpz_powm(v.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
        return;
    }
    //std::cerr<<"here2\n";
    m = p-1;
    unsigned long s = mpz_scan1(m.get_mpz_t(), 0), rk;  // s è l'esponente di 2 più grande in p-1
    //std::cerr<<s<<"\n"<<std::flush;
    m >>= s;     // adesso m è la parte dispari di p-1
    //std::cerr<<s<<" "<<m<<"\n"<<std::flush;
    mpz_powm(u.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());   // u == a^m mod p
    m = (m+1)/2;
    mpz_powm(v.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());   // v == a^((m+1)/2) mod p
    m = 2*m - 1;
    //std::cerr<<u<<" "<<v<<" " <<(v*v - a*u) %p <<"\n"<<std::flush;
    if(u==1) return;

    getNonQResidue(p, c);
    mpz_powm(c.get_mpz_t(), c.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());   // da qui m non dovrebbe servire più; c è nel sylow subgroup

    while(u != 1){
        rk = getOrdInSylow(p, u, s, m);
        m = 1;
        m <<= (s-rk-1); // m = 2^(s-rk-1)
        //mpz_pow_ui(m.get_mpz_t(), m.get_mpz_t(), s-rk-1);

        mpz_powm(c.get_mpz_t(), c.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
        v = (v*c) % p;
        c = (c*c) % p;
        u = (u*c) % p;

        s = rk; // i nuovi elementi hanno ord che divide 2^rk
    }
}



unsigned long getOrdInSylow(const mpz_class& p, mpz_class& u, unsigned long s, mpz_class& tmp){
    unsigned long rk = 0;   // esponente da dare a 2 per avere u^(2^rk) == 1 % p
    tmp = u;
    while(tmp != 1 && rk < s){
        tmp = (tmp*tmp)%p;
        rk++;
    }
    return rk;
}

void getNonQResidue(const mpz_class& p, mpz_class& c){
    //std::cerr<<"here\n"<<std::flush;
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    //std::cerr<<"here\n"<<std::flush;
    int l;
    do{
        //std::cerr<<"here "<<p<<"\t"<<std::flush;
        mpz_urandomm(c.get_mpz_t(), state, p.get_mpz_t());
        //std::cerr<<c<<"\t"<<std::flush;
        l = mpz_legendre(c.get_mpz_t(), p.get_mpz_t());
        //std::cerr<<l<<"\n"<<std::flush;

    } while(l != -1);

    gmp_randclear(state);
}