#include "print_stuff.hpp"



void printFactorBase(std::unordered_map<mpz_class, unsigned short>& F){
    std::cerr << F.size() << std::endl << std::endl;
    for(auto coppia : F){
        std::cerr << "{" << coppia.first << ", " << coppia.second << "}," << std::endl;
    }
    std::cerr << std::endl;
}


void printSetaccio(std::unordered_map<mpz_class, elemSetaccio>& S, mpz_class& base){

    std::map<mpz_class, elemSetaccio*> orderedS;

    for (auto i = S.begin(); i != S.end(); i++){
        orderedS[i->first] = &(i->second);
    }

    for (auto i = orderedS.begin(); i != orderedS.end(); i++){
        std::cerr << i->first + base << ":\n\t";
        for(auto j = (i->second)->begin(); j != (i->second)->end(); j++){
            std::cerr << "[" << j->first << ", " << (double(j->second) / (1<<6)) << "]  ";
        }
        std::cerr << "\n";
    }
    std::cerr << std::endl;
}


void printSmooths(std::vector<smoothElem> &smooths){
    for(auto& coso : smooths){
        std::cerr << coso.x << " " << coso.y << ":\n\t";
        for(auto p : coso.primes){
            std::cerr << p.first << " ";
        }
        std::cerr << std::endl;
    }

    std::cerr << "(Almost) smooth elements: " << smooths.size() << std::endl << std::endl;
}



void printUsedPrimes(std::unordered_map<mpz_class, uint32_t> & S){
    std::set<mpz_class> orderedS;

    for (auto& i : S){
        orderedS.insert(i.first);
    }

    std::cerr << "Total number of primes: " << S.size() << " = " << orderedS.size() << std::endl;
    std::cerr << "Used primes: \n";
    for (auto i : orderedS){
        std::cerr << i << " ";
    }
    std::cerr << std::endl << std::endl;;
}


void printIdeal(ideal& i){
    std::cout << "{" << i.p << ", " << i.r << ", " << int(i.logP) << "}";
}


void printFBgnfs(factorBase& FB){
    std::cout << "[" << std::endl;
    for(auto& ideale : FB){
        std::cout << "\t";
        printIdeal(ideale);
        std::cout << std::endl;
    }
    std::cout << "]" << std::endl << std::endl;
}