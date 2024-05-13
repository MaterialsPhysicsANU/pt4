#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <limits>
#include <cmath>

#include "vec.h"

#ifdef _MCS_VER
	#define OS_MAX_PATH _MAX_PATH
#else
	#define OS_MAX_PATH FILENAME_MAX
#endif


extern double pi;

enum plane{
    XY,
    YZ,
    ZX,
};

enum print_mode{
    DBG,
    INF,
    WRN,
    ERR,
};

template <typename T> 
constexpr inline T sgn(const T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> 
constexpr inline T sq(const T val) {
    return val*val;
}

std::string operator/(const std::string& a, const std::string& b);

class logger{
    public:
        print_mode level;

        std::ostream* printer;

        logger(print_mode, std::ostream& = std::cout);

        std::ostream& log(print_mode);
};

extern logger out;
extern std::ostream nullstream;

int preprocess(const char*,std::FILE*);

#endif
