#include "utils_ZZX.hpp"

void ZZX_eval(ZZ& res, const ZZX& f, const ZZ& var){
    res = 0;
    thread_local ZZ pot;
    pot = 1;                        // separato dalla dichiarazione altrimenti avverrebbe una volta sola
    for(long i=0; i<=deg(f); i++){
        res += (coeff(f, i)*pot);
        pot*=var;
    }
}