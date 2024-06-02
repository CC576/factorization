#include "print_stuff.hpp"



void printFactorBase(std::vector<std::pair<mpz_class, unsigned short>>& F){
    std::cerr << std::endl << std::endl << F.size() << std::endl << std::endl;
    for(auto coppia : F){
        std::cerr << coppia.first << std::endl;
    }
    std::cerr << std::endl;
}


void printSetaccio(std::unordered_map<mpz_class, elemSetaccio>& S, mpz_class& base){

    std::map<mpz_class, elemSetaccio*> orderedS;

    for (auto i = S.begin(); i != S.end(); i++){
        orderedS[i->first + base] = &(i->second);
    }

    for (auto i = orderedS.begin(); i != orderedS.end(); i++){
        std::cerr << i->first << ":\n\t";
        for(auto j = (i->second)->begin(); j != (i->second)->end(); j++){
            std::cerr << "[" << j->first << ", " << (double(j->second) / (1<<6)) << "]  ";
        }
        std::cerr << "\n";
    }
    std::cerr << std::endl;
}