#ifndef UTILS_ZZX_HPP
#define UTILS_ZZX_HPP

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
using namespace NTL;

void ZZX_eval(ZZ& res, const ZZX& f, const ZZ& var);

#endif