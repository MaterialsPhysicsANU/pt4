#include <iostream>

#include <chrono>
#include <filesystem>
#include <omp.h>


#include "cmd_actions.h"


std::chrono::time_point<std::chrono::high_resolution_clock> tik,tok;

int test_vol_perf_file(std::ostream& os,const std::string& pt4_name){

        std::string pt4_fp = pt4_name + ".pt4";//input file
        std::string ncf_dir = pt4_name + '/';//output directory

        pt4 pt4_0 = load_pt4(pt4_fp);
        std::filesystem::remove_all(ncf_dir);
        std::filesystem::create_directories(ncf_dir / "vol");

        std::chrono::time_point<std::chrono::high_resolution_clock> tik,tok;
        os << "STARTING " << "VOLUMES" << "..." << std::endl;
        tik = std::chrono::high_resolution_clock::now();

        write_volumes(pt4_name, pt4_0, false);

        tok = std::chrono::high_resolution_clock::now();
        os << "FINISHED " << "VOLUMES" << ". TIME: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik).count() << " ms" << std::endl;

        return 0;
}

int test_proj_perf_file(std::ostream& os,const std::string& pt4_name){

        std::string pt4_fp = pt4_name + ".pt4";//input file
        std::string ncf_dir = pt4_name + '/';//output directory

        pt4 pt4_0 = load_pt4(pt4_fp);
        std::filesystem::remove_all(ncf_dir);
        std::filesystem::create_directories(ncf_dir / "proj");

        std::chrono::time_point<std::chrono::high_resolution_clock> tik,tok;
        os << "STARTING " << "PROJECTIONS" << "..." << std::endl;
        tik = std::chrono::high_resolution_clock::now();

        write_projections(pt4_name, pt4_0);

        tok = std::chrono::high_resolution_clock::now();
        os << "FINISHED " << "PROJECTIONS" << ". TIME: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik).count() << " ms" << std::endl;

        return 0;
}

int test_vol_perf(){
    int ret;
    #if  defined(__unix__)
    const std::string iodir = "/g/data/w09/sc3563/patch_sensitivity/";
    #elif defined(WIN32)
    const std::string iodir = "../verification/sim/perf/";
    #endif

    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(iodir)) {
        if (entry.is_regular_file()) {
            files.push_back( iodir / entry.path().filename().stem().string() );
            out.log(INF) << iodir / entry.path().filename().stem().string() << std::endl;
        }
    }
    // std::reverse(files.begin(), files.end());

    auto of = std::ofstream(iodir/"overlay"/"vol_perf.log");
    for(size_t f = 0; f < files.size(); ++f){
        const std::string fp = files[f];
        std::filesystem::create_directory(fp);
        of << "running: " << fp << std::endl;
        test_vol_perf_file(of,fp);
    }
    of.close();

    return 0;
}


int test_proj_perf(){
    int ret;
    #if  defined(__unix__)
    const std::string iodir = "/g/data/w09/sc3563/patch_sensitivity/";
    #elif defined(WIN32)
    const std::string iodir = "../verification/sim/perf/";
    #endif

    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(iodir)) {
        if (entry.is_regular_file()) {
            files.push_back( iodir / entry.path().filename().stem().string() );
            out.log(INF) << iodir / entry.path().filename().stem().string() << std::endl;
        }
    }
    // std::reverse(files.begin(), files.end());

    auto of = std::ofstream(iodir/"overlay"/"proj_perf.log");
    for(size_t f = 0; f < files.size(); ++f){
        const std::string fp = files[f];
        std::filesystem::create_directory(fp);
        of << "running: " << fp << std::endl;
        test_proj_perf_file(of,fp);
    }
    of.close();

    return 0;
}

int main(int argc, char* argv[]){
#if _OPENMP
    omp_set_max_active_levels(8);
#endif

    if (!(argc > 1)) {
        std::cerr << "fatal error: no test selected" << std::endl;
        return -2;
    }
    std::string arg = argv[1];

    if(arg == "vol_perf"){ return test_vol_perf();}
    if(arg == "proj_perf"){ return test_proj_perf();}

    std::cerr << "fatal error: invalid test" << std::endl;

    return -1;
}
