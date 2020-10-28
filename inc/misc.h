#ifndef MISC_H
#define MISC_H

#include <ctime>
#include <string>
#include "basic.h"

FLM_SCALAR duration(const clock_t &startTime, const clock_t &endTime);

void runtime_str(std::string &ret);

#endif
