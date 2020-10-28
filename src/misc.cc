#include <chrono>
#include <sstream>
#include <iomanip>
#include "../inc/misc.h"

FLM_SCALAR duration(const clock_t &startTime, const clock_t &endTime)
{
    return static_cast<FLM_SCALAR>(endTime - startTime) / CLOCKS_PER_SEC;
}

void runtime_str(std::string &ret)
{
    auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream ss;
    ss << std::put_time(std::localtime(&tt), "%Y%m%d-%H%M%S");
    ret = ss.str();
}
