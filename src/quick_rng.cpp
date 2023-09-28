
#include "quick_rng.h"

#include <utility>

void xorshift::srand(uint32_t seed){
    state = seed;
}


uint32_t xorshift::rand(uint32_t shifts){
    for(; shifts > 0; shifts--){
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
    }
    return state;
}


inline void noise_volume::srand(uint32_t s){
    seed = s;
    volume_cache.clear();
}

inline std::pair<int, uint16_t> floor_mod(int d, uint16_t m) {
    int r = d % m;          //could be negative
    if (r < 0) { r += m; }  //but now its not
    int q = (d - r) / m;
    return std::make_pair(q, r);
}

float noise_volume::operator[](const std::array<int, 3> index){
    const std::array<std::pair<int, uint16_t>,3> pair_index = {
                                             floor_mod(index[0]-cache_chunk_size/2,cache_chunk_size),//subtract half the chunk size becase most of the time we generate around 0,0,0 so this saves us from 7 cache misses
                                             floor_mod(index[1]-cache_chunk_size/2,cache_chunk_size),
                                             floor_mod(index[2]-cache_chunk_size/2,cache_chunk_size)
                                             };
    const std::array<int,3> cache_index = {
                                    pair_index[0].first,
                                    pair_index[1].first,
                                    pair_index[2].first,
                                    };
    const std::array<size_t,3> chunk_index = {
                                    pair_index[0].second,
                                    pair_index[1].second,
                                    pair_index[2].second,
                                    };

    if(volume_cache.contains(cache_index)){
        return volume_cache[cache_index][chunk_index];
    }else{
        volume_cache[cache_index] = gen_chunk(cache_index);
        return volume_cache[cache_index][chunk_index];
    }
}

inline uint32_t rotr32 (uint32_t n, unsigned int c){
  const unsigned int mask = (CHAR_BIT*sizeof(n) - 1);
  c &= mask;
  return (n>>c) | (n<<( (-c)&mask ));
}

volume<float> noise_volume::gen_chunk(std::array<int,3> cache_index){
    auto chunk = volume<float>({cache_chunk_size,cache_chunk_size,cache_chunk_size});

    //lets use the seed and some bit packing to get the entropy right
    xorgen.srand(   seed
                    ^rotr32(cache_index[0],0*CHAR_BIT)
                    ^rotr32(cache_index[1],1*CHAR_BIT)
                    ^rotr32(cache_index[2],2*CHAR_BIT)
                ); 
    //run xorshift a few times to make sure we are in a good range
    xorgen.rand(37u);
    std::for_each(chunk.data.begin(), chunk.data.end(), [this](float &v){v = 2*(xorgen.rand(1u)/(float)UINT32_MAX)-1.0;});
    return chunk; 
}

thread_local noise_volume texture_noise_volume = noise_volume();

