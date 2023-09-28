#ifndef SIM_H
#define SIM_H

#include "util.h"
#include "ndvec.h"
#include "pt4.h"
#include "primative.h"
#include <memory>
#include <random>

int add_noise_projection(slice<float>& projection, const scan_spec_t& scan_spec, std::default_random_engine& rng);
int project_pt4(slice<float>& projection, const pt4& pt4_0, double time, double angle, unsigned int samples);
void add_pt4_vol(volume<float> &vol, const pt4& pt4_0,const int t_frame, const uint8_t msaa);

#endif