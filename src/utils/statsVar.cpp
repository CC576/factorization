#include "statsVar.hpp"

#include<iostream>
#include<sys/resource.h>

bool printStats = false;

int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
/* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}


void printElapsedTime(const std::string& msg, timeval& t1){   // assume che t1 contenga il precedente tempo usato
                                                              // e ci mette quello nuovo
    thread_local struct rusage usage;
    thread_local timeval t2, tdiff;

    getrusage(RUSAGE_SELF, &usage);
    t2 = usage.ru_utime;
    timeval_subtract(&tdiff, &t2, &t1);

    t1 = t2;

    std::clog << "\"" << msg << "\": [" << tdiff.tv_sec << ", " << tdiff.tv_usec << "], ";
}


void initT1(timeval& t1){
  thread_local struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  t1 = usage.ru_utime;
}