#ifndef STATS_VAR_HPP
#define STATS_VAR_HPP

#include<time.h>
#include <string>

extern bool printStats;
int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);
void initT1(timeval& t1);
void printElapsedTime(const std::string& msg, timeval& t1);

#endif