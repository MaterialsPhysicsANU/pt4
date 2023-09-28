#ifndef SCAN_SPEC_H
#define SCAN_SPEC_H

#ifdef __cplusplus
    #include "util.h"
#endif

typedef enum{
	INIT_PROJ = 0,
	INTENSITY,      //solid fill
	ATTENUATION,    //alternating black white planes
}proj_type;

#ifdef __cplusplus
    int uninitalised_proj_type();
#endif

struct scan_spec_t{
    proj_type proj_int;

    int seed;
    double I0;//photon flux (photons/projection/[0,1]^2)
    int noise_quanisation;//(bool)
    int noise_poisson;//(bool)
    double noise_gaussian;//(stdev photons/pixel)

    double revolutions_per_unit_time;
    double projections_per_revolution;
    unsigned int projection_supersampling_ratio;

};

#endif