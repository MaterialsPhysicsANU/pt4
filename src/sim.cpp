
#include "sim.h"

#include "scan_spec.h"
#include "pt4.h"

#include <random>
#include <cmath>
#include "quick_rng.h"

float blend(float a,float b, float alpha, blend_type blend){
    float c;
    switch (blend)
    {
    case REPLACE:
        c = b;
        break;
    case ADD:
        c = a+b;
        break;
    case MULTIPLY:
        c = a*b;
        break;
    case MASK:
        c = (b > 0) ? b : a;
        break;
    default:
        return NAN;
        break;
    }
    return c*alpha + a*(1-alpha);
}

//psudo-psudorandom gen
std::array<double,3> voxel_rng(uint32_t s){
    if(s == 0){return {0.5,0.5,0.5};}
    xorshift msaa_rng; msaa_rng.srand(s);
    double subx = ((msaa_rng.rand(6u)%UINT8_MAX)/(double)UINT8_MAX);
    double suby = ((msaa_rng.rand(6u)%UINT8_MAX)/(double)UINT8_MAX);
    double subz = ((msaa_rng.rand(6u)%UINT8_MAX)/(double)UINT8_MAX);

    std::array<double,3> subdisp = {subx,suby,subz};
    return subdisp;
};

int add_noise_projection(slice<float>& projection, const scan_spec_t& scan_spec, std::default_random_engine& rng){
    const double pixel_area = 4.0/product(projection.size);//4.0 because we consider projection to have dims [-1,1]^2.
    switch(scan_spec.proj_int){
            case INIT_PROJ:
                exit(uninitalised_proj_type());
                break;
            case INTENSITY:
                for(float& pv : projection.data){
                    const double photons = pixel_area*scan_spec.I0 * exp(-((double)pv/pixel_area));//by doing this division we always exp near the actual intensity level
                    double sample_photons = 0.0;
                    if(!(scan_spec.noise_gaussian == 0.0)){
                        std::normal_distribution<double> gaussian(0.0, scan_spec.noise_gaussian);
                        sample_photons += gaussian(rng);
                    }
                    if(scan_spec.noise_poisson){
                        std::poisson_distribution<int> poisson(photons);
                        sample_photons += poisson(rng);
                    }else{
                        sample_photons += photons;
                    }
                    if(scan_spec.noise_quanisation){
                        sample_photons = std::round(sample_photons);
                    }
                    pv = sample_photons;
                }
                break;
            case ATTENUATION:
                break;
    }
    return 0;
}


int project_pt4(slice<float>& projection, const pt4& pt4_0, double time, double angle, unsigned int samples){
    #pragma omp parallel
	{
	location camera;
    camera.angle = angle; 

    const float voxel_volumei = 8.0/(projection.size[0]*projection.size[1]*samples);//8 because [-1, 1]^3
    auto primatives = std::vector<std::unique_ptr<primative>>();
    for (lazy_primative l_p : pt4_0.primitives) {
        if (l_p.domains_h->t_length() >= time) {
            std::unique_ptr<primative> p = primative_from_lazy(l_p, time);
            primatives.push_back(std::move(p));
        }
    }

	#pragma omp for collapse(2)
    for(size_t it = 0; it < projection.size[0]; it ++){
    for(size_t iz = 0; iz < projection.size[1]; iz ++){
        auto vox = std::array<double, 2> {(double)it+0.5,(double) iz+0.5};

        std::array<double,2> detector = voxel2screen(projection.size,vox);//convert to screen space


        float integral = 0.0;
        for(unsigned int s = 0; s < samples; s++){
            std::array<double,3> screen = {detector[0], (2.0*(s+0.5))/samples-1.0, detector[1]};//+0.5 to get center of voxel

            screen = camera.inv_affine(screen);

            float a = 0.0;
            float b;
            for(std::unique_ptr<primative>& p : primatives){
                std::array<double,3> prim_space = p->loc.inv_affine(screen);
                if(p->inside(prim_space)){
                    b = p->atn.atn_at(prim_space);
                    a = blend(a,b,1,p->atn.blend);
                }
            }
            integral += a;
        }
        projection[{it,iz}] = integral*voxel_volumei;
    }
    }

	}
    return 0;
}

void add_pt4_vol(volume<float> &vol, const pt4& pt4_0, int t_frame, const uint8_t msaa){
    #pragma omp parallel
	{
    const float voxel_volumei = 8.0/(vol.size[0]*vol.size[1]*vol.size[2]);//8 because [-1, 1]^3
    auto primatives = std::vector<std::unique_ptr<primative>>();
    for (lazy_primative l_p : pt4_0.primitives) {
        double time = t_frame*pt4_0.vol_t_step;
        if (l_p.domains_h->t_length() >= time) {
            std::unique_ptr<primative> p = primative_from_lazy(l_p, time);
            primatives.push_back(std::move(p));
        }
    }

	#pragma omp for collapse(3)
	for(size_t ivox = 0; ivox < vol.size[0]; ivox++){
    for(size_t jvox = 0; jvox < vol.size[1]; jvox++){
    for(size_t kvox = 0; kvox < vol.size[2]; kvox++){
        auto vox = std::array<double, 3> {(double)ivox,(double)jvox,(double) kvox};
        
        float sample_atn_accumulant = 0.0;
        for(uint8_t s = 0; s < msaa; s++){
            std::array<double,3> vox_rng = vox + voxel_rng(s);//add msaa in voxel space
            std::array<double,3> screen = voxel2screen(vol.size,vox_rng);//convert to screen space

            float a = 0.0;
            float b;
            for(std::unique_ptr<primative>& p : primatives){
                std::array<double,3> prim_space = p->loc.inv_affine(screen);
                if(p->inside(prim_space)){
                    b = p->atn.atn_at(prim_space);
                    a = blend(a,b,1,p->atn.blend);
                }
            }
            sample_atn_accumulant += a;
        }
        vol[{ivox,jvox, kvox}] = (sample_atn_accumulant*voxel_volumei)/msaa;
        
    }
    }
    }
	
	}
}
