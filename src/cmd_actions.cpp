#define _CRT_SECURE_NO_WARNINGS

#include "cmd_actions.h"

#include <string>
#include <vector>
#include <filesystem>
#include <sstream>

#include <random>

#include "ndvec.h"
#include "pt4.h"
#include "sim.h"

#undef _CRT_SECURE_NO_WARNINGS

pt4 load_pt4(std::string& pt4_fp){
    #pragma omp critical (load_pt4)
    {
    int ret;
    std::FILE* pt4_temp_file = std::tmpfile();

    out.log(INF) << "PREPROCESSING PHANTOM..." << std::endl;

    ret = preprocess(pt4_fp.c_str(), pt4_temp_file);
    if(ret){exit(ret);}

    out.log(INF) << "PARSING PHANTOM..." << std::endl;

    auto pt4_0 = pt4(pt4_temp_file);
    fclose(pt4_temp_file);
    return pt4_0;
    }
}


int write_volumes(const std::string& odir, const pt4 pt4_0, bool zproject){
    double total_time = pt4_0.t_length();
    int total_frames = (int) (total_time/pt4_0.vol_t_step);

    auto project_x = volume<float>({pt4_0.size[1],pt4_0.size[2],(size_t)total_frames});
    auto project_y = volume<float>({pt4_0.size[2],pt4_0.size[0],(size_t)total_frames});
    auto project_z = volume<float>({pt4_0.size[0],pt4_0.size[1],(size_t)total_frames});

    int completed_frames = 0;  
    #pragma omp parallel for
    for(int time_index = 0; time_index < total_frames; time_index++){
        #pragma omp critical (wtv0)
        {
        char progress_str[0xFF]; 
	snprintf(progress_str, sizeof(progress_str), "frame: %05d/%05d\n", completed_frames, total_frames);
        out.log(INF) << progress_str;
        completed_frames++;
        }
        auto phantom = volume<float>({pt4_0.size[0],pt4_0.size[1],pt4_0.size[2]});
        add_pt4_vol(phantom, pt4_0, time_index, 8);
        
        if(zproject){
            for(size_t x = 0; x < pt4_0.size[0]; x++){
            for(size_t y = 0; y < pt4_0.size[1]; y++){
            for(size_t z = 0; z < pt4_0.size[2]; z++){
                project_x[{y,z,(size_t)time_index}] += phantom[{x,y,z}];
                project_y[{z,x,(size_t)time_index}] += phantom[{x,y,z}];
                project_z[{x,y,(size_t)time_index}] += phantom[{x,y,z}];
            }
            }
            }
        }

        char fp[OS_MAX_PATH];
        sprintf(fp, "%s/vol/vol%05d.nc", odir.c_str() ,time_index);
        vol2ncf(fp, phantom);


    }
    if(zproject){
        std::string pt4_fp;
        pt4_fp = odir + "/ztime/project_x.nc";    vol2ncf(pt4_fp.c_str(), project_x);
        pt4_fp = odir + "/ztime/project_y.nc";    vol2ncf(pt4_fp.c_str(), project_y);
        pt4_fp = odir + "/ztime/project_z.nc";    vol2ncf(pt4_fp.c_str(), project_z);
    }
    return 0;
}

int write_projections(const std::string& pt4_name, const pt4 pt4_0){
    const double total_time = pt4_0.t_length();
    const double projections_per_unit_time = pt4_0.scan_spec.revolutions_per_unit_time*pt4_0.scan_spec.projections_per_revolution;
    const int frames = total_time* projections_per_unit_time;
    const double super_samp = pt4_0.scan_spec.projection_supersampling_ratio; 

    int completed_frames = 0;
    #pragma omp parallel for
    for(int frame = 0; frame < frames; frame++){
        #pragma omp critical (wtp0)
        {
        char progress_str[0xFF]; 
        snprintf(progress_str, sizeof(progress_str), "frame: %05d/%05d\n", completed_frames, frames);
        out.log(INF) << progress_str;
        completed_frames++;
        }

        if(pt4_0.size[0] != pt4_0.size[1]){out.log(ERR) << "different x,y dims not yet supported!\n"; exit(2);}
        auto projection = slice<float>({ceil(pt4_0.size[0]*super_samp),ceil(pt4_0.size[2]*super_samp)});

        const double time  = frame/projections_per_unit_time;
        const double angle = 2.0f*pi*frame/pt4_0.scan_spec.projections_per_revolution;
        const unsigned int samples  = pt4_0.size[1]*super_samp;
        project_pt4(projection, pt4_0, time, angle, samples);

        slice<float> subsampled;
        subsampled = projection.subsample(super_samp);

        uint64_t angle_bits = std::bit_cast<uint64_t>(angle);
        uint32_t seed = (uint32_t)(angle_bits >> 0x20) ^ ((uint32_t)angle_bits) ^ pt4_0.scan_spec.seed;
        auto rng = std::default_random_engine(seed);
        add_noise_projection(subsampled, pt4_0.scan_spec, rng);


	    char fp[OS_MAX_PATH];
        sprintf(fp, "%s/proj/%05d.tiff", pt4_name.c_str(), frame);
        slice2tiff(fp, subsampled);
    }
    return 0;
}