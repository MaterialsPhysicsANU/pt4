#ifndef XORSHIFT_H
#define XORSHIFT_H

#include <cstdint>
#include <map>
#include <thread>
#include "vec.h"
#include "ndvec.h"

class xorshift {
        uint32_t state;

    public:
        void srand(uint32_t seed);
        uint32_t rand(uint32_t shifts = 1);
};

class noise_volume {
        const uint16_t cache_chunk_size = 64;
        std::map<std::array<int,3>, volume<float>> volume_cache;
        xorshift xorgen;
        uint32_t seed = 0x7D'8C'AA'FA ;

        volume<float> gen_chunk(std::array<int,3>);
    public:
        float operator[](const std::array<int, 3> index);

        inline void srand(uint32_t);
    
};

extern thread_local noise_volume texture_noise_volume;

#endif
